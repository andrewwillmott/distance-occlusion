//
// DOFImage.cpp
//
// Distance/occlusion field image IO
//
// Andrew Willmott
//

#define _CRT_SECURE_NO_WARNINGS  // for stb_image

#include "DOFImage.hpp"

#include "DistanceFields.hpp"
#include "OcclusionFields.hpp"

#include "stb_image_mini.h"


// Masks

uint32_t* DOL::ReadMask(const char* filename, int threshold, bool lessThan, int& w, int& h)
{
    stbi_uc* data = stbi_load(filename, &w, &h, 0, 1);
    if (!data)
    {
        printf("error reading %s: %s\n", filename, stbi_failure_reason());
        return 0;
    }

    const int sw = MaskSize(w);
    uint32_t* mask = new uint32_t[sw * h];

    for (int y = 0; y < h; y++)
        for (int i = 0; i < sw; i++)
        {
            uint32_t* maskCell = mask + y * sw + i;
            *maskCell = 0;

            for (int j = 0; j < 32; j++)
            {
                int x = i * 32 + j;

                if (x >= w)
                    break;

                uint8_t v = data[y * w + x];

                if ((v >= threshold) ^ lessThan)
                    *maskCell |= 1 << j;
            }
    }

    stbi_image_free(data);

    return mask;
}

namespace
{
    void WriteImage(const char* name, const char* type, int w, int h, uint8_t data[], int index = -1)
    {
        char nameAndExt[128];

        if (index >= 0)
            snprintf(nameAndExt, 128, "%s-%s-%03d.png", name, type, index);
        else
            snprintf(nameAndExt, 128, "%s-%s.png", name, type);

        stbi_write_png(nameAndExt, w, h, 1, data, 0);
        printf("written %s\n", nameAndExt);
    }

    void WriteImage(const char* name, const char* type, int w, int h, uint8_t data[][4], int index = -1)
    {
        char nameAndExt[128];

        if (index >= 0)
            snprintf(nameAndExt, 128, "%s-%s-%03d.png", name, type, index);
        else
            snprintf(nameAndExt, 128, "%s-%s.png", name, type);

        stbi_write_png(nameAndExt, w, h, 4, data, 0);
        printf("written %s\n", nameAndExt);
    }
}

void DOL::WriteMask(const char* name, int w, int h, const uint32_t mask[])
{
    uint8_t* dist8 = new uint8_t[w * h];

    for (int y = 0; y < h; y++)
    {
        uint8_t* dist8Row = dist8 + w * y;

        for (int j = 0; j < w; j += 32)
        {
            uint32_t m = *mask++;

            int n = 32;
            if (w - j < n)
                n = w - j;

            for (int i = 0; i < n; i++)
            {
                dist8Row[j + i] = (m & 1) ? 255 : 0;

                m >>= 1;
            }
        }
    }

    WriteImage(name, "mask", w, h, dist8);

    delete[] dist8;
}

// Distances

void DOL::WriteDistances(const char* name, const char* type, int w, int h, int d, float maxLen, const CellDelta3 deltas[])
{
    // Write as stack of images
    uint8_t* dist8 = new uint8_t[w * h];

    for (int slice = 0; slice < d; slice++)
    {
        const CellDelta3* sliceDeltas = deltas + w * h * slice;

        for (int i = 0, n = w * h; i < n; i++)
            dist8[i] = EncodeU8(Distance(sliceDeltas[i]) / maxLen);

        WriteImage(name, type, w, h, dist8, slice);
    }

    delete[] dist8;
}

// delta write

void DOL::WriteDistances(const char* name, int w, int h, float maxLen, const uint32_t mask[], const CellDelta2 deltas[], int imageType)
{
    if (imageType == kImageStandard)
    {
        uint8_t* dist8 = new uint8_t[w * h];

        for (int y = 0; y < h; y++)
        {
            const CellDelta2* deltasRow = deltas + w * y;
            uint8_t*          dist8Row  = dist8  + w * y;

            for (int j = 0; j < w; j += 32)
            {
                uint32_t m = *mask++;

                int n = 32;
                if (w - j < n)
                    n = w - j;

                for (int i = 0; i < n; i++)
                {
                    float d = Distance(deltasRow[j + i]) / maxLen;
                    dist8Row[j + i] = EncodeU8Signed((m & 1) ? -d : d);

                    m >>= 1;
                }
            }
        }

        WriteImage(name, "sdf", w, h, dist8);

        delete[] dist8;
    }
    else
    {
        uint8_t (*dist32)[4] = new uint8_t[w * h][4];

        for (int y = 0; y < h; y++)
        {
            const CellDelta2* deltasRow     = deltas + w * y;
            uint8_t         (*dist32Row)[4] = dist32 + w * y;

            if (imageType == kImageDeltas)
                for (int x = 0; x < w; x++)
                {
                    dist32Row[x][0] = EncodeU8Signed(deltasRow[x].x / maxLen);
                    dist32Row[x][1] = EncodeU8Signed(deltasRow[x].y / maxLen);
                    dist32Row[x][2] = 0;
                    dist32Row[x][3] = 255;
                }
            else if (imageType == kImageRG)
                for (int j = 0; j < w; j += 32)
                {
                    uint32_t m = *mask++;

                    int n = 32;
                    if (w - j < n)
                        n = w - j;

                    for (int i = 0; i < n; i++)
                    {
                        dist32Row[j + i][ m & 1     ] = EncodeU8(Distance(deltasRow[j + i]) / maxLen);
                        dist32Row[j + i][(m & 1) ^ 1] = 0;
                        dist32Row[j + i][2] = 0;
                        dist32Row[j + i][3] = 255;

                        m >>= 1;
                    }
                }
            else if (imageType == kImageAngles)
                for (int x = 0; x < w; x++)
                {
                    float d = Distance(deltasRow[x]);
                    dist32Row[x][0] = EncodeU8Signed(deltasRow[x].x / d);
                    dist32Row[x][1] = EncodeU8Signed(deltasRow[x].y / d);
                    dist32Row[x][2] = 0;
                    dist32Row[x][3] = 255;
                }
        }

        const char* suffix;
        switch (imageType)
        {
        case kImageDeltas: suffix = "sdf-delta" ; break;
        case kImageRG    : suffix = "sdf-rg"    ; break;
        case kImageAngles: suffix = "sdf-angles"; break;
        default          : suffix = "unknown"   ; break;
        }

        WriteImage(name, suffix, w, h, dist32);
        delete[] dist32;
    }
}

void DOL::WriteDistances(const char* name, int w, int h, float maxLen, const CellDelta2 deltas0[], const CellDelta2 deltas1[], int imageType)
{
    if (imageType == kImageStandard)
    {
        uint8_t* dist8 = new uint8_t[w * h];

        for (int i = 0, n = w * h; i < n; i++)
            dist8[i] = EncodeU8Signed((Distance(deltas0[i]) - Distance(deltas1[i])) / maxLen);

        WriteImage(name, "sdf", w, h, dist8);
        delete[] dist8;
    }
    else
    {
        uint8_t (*dist32)[4] = new uint8_t[w * h][4];

        for (int i = 0, n = w * h; i < n; i++)
        {
            if (imageType == kImageDeltas)
            {
                dist32[i][0] = EncodeU8Signed((deltas0[i].x + deltas1[i].x) / maxLen);
                dist32[i][1] = EncodeU8Signed((deltas0[i].y + deltas1[i].y) / maxLen);
            }
            else if (imageType == kImageRG)
            {
                dist32[i][0] = EncodeU8(Distance(deltas0[i]) / maxLen);
                dist32[i][1] = EncodeU8(Distance(deltas1[i]) / maxLen);
            }
            else if (imageType == kImageAngles)
            {
                float d = (Distance(deltas0[i]) - Distance(deltas1[i]));

                dist32[i][0] = EncodeU8Signed((deltas0[i].x + deltas1[i].x) / d);
                dist32[i][1] = EncodeU8Signed((deltas0[i].y + deltas1[i].y) / d);
            }
            dist32[i][2] = 0;
            dist32[i][3] = 255;
        }

        const char* suffix;
        switch (imageType)
        {
        case kImageDeltas: suffix = "sdf-delta" ; break;
        case kImageRG    : suffix = "sdf-rg"    ; break;
        case kImageAngles: suffix = "sdf-angles"; break;
        default          : suffix = "unknown"   ; break;
        }

        WriteImage(name, suffix, w, h, dist32);
        delete[] dist32;
    }
}

void DOL::WriteDistances(const char* name, const char* type, int w, int h, float maxLen, const CellDelta2 deltas[])
{
    uint8_t* dist8 = new uint8_t[w * h];

    for (int i = 0, n = w * h; i < n; i++)
        dist8[i] = EncodeU8(Distance(deltas[i]) / maxLen);

    WriteImage(name, type, w, h, dist8);
    delete[] dist8;
}

void DOL::WriteDistances(const char* name, const char* type, int w, int h, float maxLen, const float distances[])
{
    uint8_t* dist8 = new uint8_t[w * h];

    for (int i = 0, n = w * h; i < n; i++)
        dist8[i] = EncodeU8(distances[i] / maxLen);

    WriteImage(name, type, w, h, dist8);
    delete[] dist8;
}

void DOL::WriteDistances(const char* name, const char* type, int w, int h, float maxLen, const int16_t distances[])
{
    uint8_t* dist8 = new uint8_t[w * h];

    for (int i = 0, n = w * h; i < n; i++)
        dist8[i] = EncodeU8(distances[i] / maxLen);

    WriteImage(name, type, w, h, dist8);
    delete[] dist8;
}

void DOL::WriteDistances(const char* name, int w, int h, int d, float maxLen, const uint32_t mask[], const CellDelta3 deltas[], int imageType)
{
    if (imageType == kImageStandard)
    {
        // Write as images
        uint8_t* dist8 = new uint8_t[w * h];

        for (int slice = 0; slice < d; slice++)
        {
            for (int y = 0; y < h; y++)
            {
                const CellDelta3* deltasRow = deltas + w * y;
                uint8_t*            dist8Row = dist8  + w * y;

                for (int j = 0; j < w; j += 32)
                {
                    uint32_t m = *mask++;

                    int n = 32;
                    if (w - j < n)
                        n = w - j;

                    for (int i = 0; i < n; i++)
                    {
                        float d = Distance(deltasRow[j + i]) / maxLen;
                        dist8Row[j + i] = EncodeU8Signed((m & 1) ? -d : d);

                        m >>= 1;
                    }
                }
            }

            WriteImage(name, "sdf", w, h, dist8, slice);
            deltas += w * h;
        }

        delete[] dist8;
    }
    else
    {
        uint8_t (*dist32)[4] = new uint8_t[w * h][4];

        const char* suffix;
        switch (imageType)
        {
        case kImageDeltas: suffix = "sdf-delta" ; break;
        case kImageRG    : suffix = "sdf-rg"    ; break;
        case kImageAngles: suffix = "sdf-angles"; break;
        default          : suffix = "unknown"   ; break;
        }

        for (int slice = 0; slice < d; slice++)
        {
            for (int y = 0; y < h; y++)
            {
                const CellDelta3* deltasRow     = deltas + w * y;
                uint8_t          (*dist32Row)[4] = dist32 + w * y;

                if (imageType == kImageDeltas)
                    for (int x = 0; x < w; x++)
                    {
                        dist32Row[x][0] = EncodeU8Signed(deltasRow[x].x / maxLen);
                        dist32Row[x][1] = EncodeU8Signed(deltasRow[x].y / maxLen);
                        dist32Row[x][2] = EncodeU8Signed(deltasRow[x].z / maxLen);
                        dist32Row[x][3] = 255;
                    }
                else if (imageType == kImageRG)
                    for (int j = 0; j < w; j += 32)
                    {
                        uint32_t m = *mask++;

                        int n = 32;
                        if (w - j < n)
                            n = w - j;

                        for (int i = 0; i < n; i++)
                        {
                            dist32Row[j + i][ m & 1     ] = EncodeU8(Distance(deltasRow[j + i]) / maxLen);
                            dist32Row[j + i][(m & 1) ^ 1] = 0;
                            dist32Row[j + i][2] = 0;
                            dist32Row[j + i][3] = 255;

                            m >>= 1;
                        }
                    }
                else if (imageType == kImageAngles)
                    for (int x = 0; x < w; x++)
                    {
                        float d = Distance(deltasRow[x]);
                        dist32Row[x][0] = EncodeU8Signed(deltasRow[x].x / d);
                        dist32Row[x][1] = EncodeU8Signed(deltasRow[x].y / d);
                        dist32Row[x][2] = EncodeU8Signed(deltasRow[x].z / d);
                        dist32Row[x][3] = 255;
                    }
            }

            WriteImage(name, suffix, w, h, dist32, slice);
            deltas += w * h;
        }

        delete[] dist32;
    }
}

void DOL::WriteDistances(const char* name, int w, int h, int d, float maxLen, const CellDelta3 deltas0[], const CellDelta3 deltas1[], int imageType)
{
    if (imageType == kImageStandard)
    {
        uint8_t* dist8 = new uint8_t[w * h];

        for (int slice = 0; slice < d; slice++)
        {
            const CellDelta3* sliceDeltas0 = deltas0 + w * h * slice;
            const CellDelta3* sliceDeltas1 = deltas1 + w * h * slice;

            for (int i = 0, n = w * h; i < n; i++)
                dist8[i] = EncodeU8Signed((Distance(sliceDeltas0[i]) - Distance(sliceDeltas1[i])) / maxLen);

            WriteImage(name, "sdf", w, h, dist8, slice);
        }

        delete[] dist8;
    }
    else
    {
        uint8_t (*dist32)[4] = new uint8_t[w * h][4];

        const char* suffix;
        switch (imageType)
        {
        case kImageDeltas: suffix = "sdf-delta" ; break;
        case kImageRG    : suffix = "sdf-rg"    ; break;
        case kImageAngles: suffix = "sdf-angles"; break;
        default          : suffix = "unknown"   ; break;
        }

        for (int slice = 0; slice < d; slice++)
        {
            const CellDelta3* sliceDeltas0 = deltas0 + w * h * slice;
            const CellDelta3* sliceDeltas1 = deltas1 + w * h * slice;

            if (imageType == kImageDeltas)
                for (int i = 0, n = w * h; i < n; i++)
                {
                    dist32[i][0] = EncodeU8Signed((sliceDeltas0[i].x + sliceDeltas1[i].x) / maxLen);
                    dist32[i][1] = EncodeU8Signed((sliceDeltas0[i].y + sliceDeltas1[i].y) / maxLen);
                    dist32[i][2] = EncodeU8Signed((sliceDeltas0[i].z + sliceDeltas1[i].z) / maxLen);
                    dist32[i][3] = 255;
                }
            else if (imageType == kImageRG)
                for (int i = 0, n = w * h; i < n; i++)
                {
                    dist32[i][0] = EncodeU8(Distance(sliceDeltas0[i]) / maxLen);
                    dist32[i][1] = EncodeU8(Distance(sliceDeltas1[i]) / maxLen);
                    dist32[i][2] = 0;
                    dist32[i][3] = 255;
                }
            else if (imageType == kImageAngles)
                for (int i = 0, n = w * h; i < n; i++)
                {
                    float d = (Distance(deltas0[i]) - Distance(deltas1[i]));

                    dist32[i][0] = EncodeU8Signed((deltas0[i].x + deltas1[i].x) / d);
                    dist32[i][1] = EncodeU8Signed((deltas0[i].y + deltas1[i].y) / d);
                    dist32[i][2] = EncodeU8Signed((deltas0[i].z + deltas1[i].z) / d);
                    dist32[i][3] = 255;
                }

            WriteImage(name, suffix, w, h, dist32, slice);
        }

        delete[] dist32;
    }
}

void DOL::WriteOcclusion(const char* name, int w, int h, const float dirW[])
{
    const size_t imageSize = w * h;
    const float* dirW4[4] = { dirW, dirW + imageSize, dirW + imageSize * 2, dirW + imageSize * 3 };

    uint8_t* ao = new uint8_t[w * h];

    for (int i = 0, n = w * h; i < n; i++)
    {
        const float dirWCell[4] = { dirW4[0][i], dirW4[1][i], dirW4[2][i], dirW4[3][i] };

        ao[i] = EncodeAO4(dirWCell);
    }

    WriteImage(name, "ao", w, h, ao);

    uint8_t (*aoDir)[4] = new uint8_t[w * h][4];

    for (int i = 0, n = w * h; i < n; i++)
    {
        const float dirWCell[4] = { dirW4[0][i], dirW4[1][i], dirW4[2][i], dirW4[3][i] };

        EncodeDirW4(dirWCell, aoDir[i]);
   }

    WriteImage(name, "dirW", w, h, aoDir);

    for (int d = 0; d < 4; d++)
    {
        for (int i = 0, n = w * h; i < n; i++)
            ao[i] = EncodeU8(1.0f - dirW4[d][i]);

        WriteImage(name, "dirW", w, h, ao, d);
    }

    delete[] ao;
    delete[] aoDir;
}

void DOL::WriteOcclusion(const char *name, int w, int h, int d, const float dirW[])
{
    // Write as images
    const size_t volumeStride = w * h * d;
    const size_t sliceStride = w * h;

    uint8_t* ao         = new uint8_t[sliceStride];
    uint8_t (*aoDir)[4] = new uint8_t[sliceStride][4];

    for (int slice = 0; slice < d; slice++)
    {
        const float* sDirW[8];

        for (int d = 0; d < 8; d++)
            sDirW[d] = dirW + sliceStride * slice + volumeStride * d;

        for (int i = 0, n = w * h; i < n; i++)
        {
            const float dirWCell[8] = { sDirW[0][i], sDirW[1][i], sDirW[2][i], sDirW[3][i],
                                        sDirW[4][i], sDirW[5][i], sDirW[6][i], sDirW[7][i] };

            ao[i] = EncodeAO8(dirWCell);
        }

        WriteImage(name, "ao", w, h, ao, slice);

        for (int i = 0, n = w * h; i < n; i++)
        {
            const float dirWCell[8] = { sDirW[0][i], sDirW[1][i], sDirW[2][i], sDirW[3][i],
                                        sDirW[4][i], sDirW[5][i], sDirW[6][i], sDirW[7][i] };

            EncodeDirW8(dirWCell, aoDir[i]);
        }

        WriteImage(name, "dirW", w, h, aoDir, slice);
    }

    delete[] ao;
    delete[] aoDir;
}

void DOL::WriteShadow(const char* name, const char* type, int w, int h, const float shadow[], int index)
{
    int n = w * h;
    uint8_t* occ8 = new uint8_t[n];

    for (int i = 0; i < n; i++)
        occ8[i] = EncodeU8(1.0f - shadow[i]);

    WriteImage(name, "shadow", w, h, occ8, index);
    delete[] occ8;
}