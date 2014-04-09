//------------------------------------------------------------------------------
// Purpose: Distance Field and Occlusion Lib test app
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#include "DistanceOcclusionLib.h"

#include <stdio.h>
#include <math.h>
#include <float.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifndef USE_TGA
    #define USE_TGA 1
#endif

#ifndef USE_OBJ
    #define USE_OBJ 1
#endif

#if USE_TGA
    #include "targa.h"
#endif

#if USE_OBJ
    #include "MeshSupport.h"
    using namespace MSL;
#endif

using namespace DOL;

namespace
{
    void PrintMask(int w, int h, const uint32_t* mask, FILE* s)
    {
        int maskStride = (w + 31) / 32;
        
        for (int y = 0; y < h; y++)
        {
            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m = *mask++;
                            
                int n = 32;
                if (w - j < n)
                    n = w - j;
                        
                for (int i = 0; i < n; i++)
                {
                    if (m & 1)
                        fprintf(s, "X");
                    else
                        fprintf(s, ".");
                    
                    m >>= 1;
                }
            }
            
            fprintf(s, "\n");
        }
    }

    void PrintMask(int w, int h, int d, const uint32_t mask[], FILE* s)
    {
        const int maskStride = (w + 31) / 32;
        
        for (int z = 0; z < d; z++)
        {
            for (int y = 0; y < h; y++)
            {
                for (int j = 0; j < maskStride * 32; j += 32)
                {
                    uint32_t m = (*mask++);
                    
                    int n = 32;
                    if (w - j < n)
                        n = w - j;
                    
                    for (int i = 0; i < n; i++)
                    {
                        if (m & 1)
                            fprintf(s, "X");
                        else
                            fprintf(s, ".");

                        m >>= 1;
                    }
                }
                
                fprintf(s, "\n");
            }
            
            fprintf(s, "\n");
        }
    }


    void PrintDistances(int w, int h, const float* dist, FILE* s)
    {    
        for (int i = 0; i < h; i++)
        {
            for (int i = 0; i < w - 1; i++)
                fprintf(s, "%1.3f, ", (*dist++));
            
            fprintf(s, "%1.3f\n", (*dist++));
        }
    }

    void PrintDistances(int w, int h, const int32_t* dist, FILE* s)
    {    
        for (int i = 0; i < h; i++)
        {
            for (int i = 0; i < w - 1; i++)
                fprintf(s, "%4d, ", (*dist++));
            
            fprintf(s, "%4d\n", (*dist++));
        }
    }

    void PrintDistances(int w, int h, int d, const float* dist, FILE* s)
    {    
        for (int k = 0; k < d; k++)
        {
            fprintf(s, "slice %d>\n", k);
            
            for (int j = 0; j < h; j++)
            {
                for (int i = 0; i < w - 1; i++)
                    fprintf(s, "%1.3f, ", (*dist++));
                    
                fprintf(s, "%1.3f\n", (*dist++));
            }
        }
    }

    void PrintDistances(int w, int h, int d, const int32_t* dist, FILE* s)
    {    
        for (int k = 0; k < d; k++)
        {
            fprintf(s, "slice %d>\n", k);
            
            for (int j = 0; j < h; j++)
            {
                for (int i = 0; i < w - 1; i++)
                    fprintf(s, "%4d, ", (*dist++));
                
                fprintf(s, "%4d\n", (*dist++));
            }
        }
    }

    void PrintDeltas(int w, int h, const cCellDelta2* deltas, FILE* s)
    {    
        for (int i = 0; i < h; i++)
        {
            for (int i = 0; i < w - 1; i++, deltas++)
                fprintf(s, "%3d:%3d, ", deltas->x, deltas->y);
            
            fprintf(s, "%3d:%3d\n", deltas->x, deltas->y);
            deltas++;
        }
    }

    void PrintDeltas(int w, int h, int d, const cCellDelta3* deltas, FILE* s)
    {    
        for (int k = 0; k < d; k++)
        {
            fprintf(s, "slice %d>\n", k);
            
            for (int j = 0; j < h; j++)
            {
                for (int i = 0; i < w; i++, deltas++)
                {
                    fprintf(s, "%3d %3d %3d", deltas->x, deltas->y, deltas->z);
                
                    if (i < w - 1)
                        fprintf(s, "|");
                    else
                        fprintf(s, "\n");
                }
            }
        }
    }

    void PrintNonZeroDeltas(int w, int h, const cCellDelta2* deltas, FILE* s)
    {
        for (int i = 0; i < h; i++)
        {
            for (int j = 0; j < w; j++, deltas++)
            {
                if (deltas->x == 0 && deltas->y == 0)
                    fprintf(s, "   *   ");
                else
                    fprintf(s, "%3d:%3d", deltas->x, deltas->y);
                
                if (j < w - 1)
                    fprintf(s, ", ");
                else
                    fprintf(s, "\n");
            }
        }
    }

    void PrintDiff(int w, int h, const cCellDelta2* deltas, const cCellDelta2* refs, FILE* s)
    {
        float maxDist = 0.0f;

        for (int i = 0; i < h; i++)
        {
            for (int j = 0; j < w; j++, deltas++, refs++)
            {
                if (deltas->x == refs->x && deltas->y == refs->y)
                    fprintf(s, "   -   ");
                else
                {
                    int dd2 = deltas->x * deltas->x + deltas->y * deltas->y;
                    int rd2 = refs  ->x * refs  ->x + refs  ->y * refs  ->y;
                    
                    float dist = sqrtf(dd2) - sqrtf(rd2);
                    
                    if (dist < 1e-6f)
                        fprintf(s, "   x   ");
                    else
                    {
                        fprintf(s, "%7.2f", dist);
                        if (maxDist < dist)
                            maxDist = dist;
                    }
                }
                
                if (j < w - 1)
                    fprintf(s, ", ");
                else
                    fprintf(s, "\n");
            }
        }

        fprintf(s, "Max error: %g\n", maxDist);
        if (maxDist >= 1.0f)
            fprintf(s, "FAILED\n");
    }

    void PrintDiff(int w, int h, int d, const cCellDelta3* deltas, const cCellDelta3* refs, FILE* s)
    {
        float maxDist = 0.0f;
        
        for (int k = 0; k < d; k++)
        {
            fprintf(s, "slice %d>\n", k);
            
            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++, deltas++, refs++)
                {
                    if (deltas->x == refs->x && deltas->y == refs->y && deltas->z == refs->z)
                        fprintf(s, "   -   ");
                    else
                    {
                        int dd2 = deltas->x * deltas->x + deltas->y * deltas->y + deltas->z * deltas->z;
                        int rd2 = refs  ->x * refs  ->x + refs  ->y * refs  ->y + refs  ->z * refs  ->z;
                        
                        float dist = sqrtf(dd2) - sqrtf(rd2);
                        
                        if (dist < 1e-6f)
                            fprintf(s, "   x   ");
                        else
                        {
                            fprintf(s, "%7.2f", dist);
                            
                            if (maxDist < dist)
                                maxDist = dist;
                        }
                    }
                    
                    if (j < w - 1)
                        fprintf(s, ", ");
                    else
                        fprintf(s, "\n");
                }
            }
        }
        
        fprintf(s, "Max error: %g\n", maxDist);
        if (maxDist >= 1.0f)
            fprintf(s, "FAILED\n");
    }

    bool CheckDeltas(int w, int h, cCellDelta2* delta, FILE* s)
    {
        bool result = true;
        
        for (int j = 0; j < h; j++)
            for (int i = 0; i < w; i++)
            {
                cCellDelta2 cell = delta[i + w * j];
                cell.x += i;
                cell.y += j;
                
                if (cell.x < 0 || cell.x >= w || cell.y < 0 || cell.y >= h)
                {
                    fprintf(s, "out of bounds cp (%d:%d) @ cell %d:%d\n", cell.x, cell.y, i, j);
                    result = false;
                    continue;
                }
                
                cCellDelta2 scell = delta[cell.x + w * cell.y];
                
                if (scell.x != 0 || scell.y != 0)
                {
                    result = false;
                    fprintf(s, "non-boundary cp (%d:%d) @ cell %d:%d\n", scell.x, scell.y, i, j);
                }
            }
            
        return result;
    }

    bool CheckDeltas(int w, int h, int d, cCellDelta3* delta, FILE* s)
    {
        bool result = true;
        
        for (int k = 0; k < d; k++)
            for (int j = 0; j < h; j++)
                for (int i = 0; i < w; i++)
                {
                    cCellDelta3 cell = delta[i + w * j + w * h * k];
                    
                    cell.x += i;
                    cell.y += j;
                    cell.z += k;
                    
                    if (cell.x < 0 || cell.x >= w || cell.y < 0 || cell.y >= h || cell.z < 0 || cell.z >= d)
                    {
                        fprintf(s, "out of bounds cp (%d:%d:%d) @ cell %d:%d:%d\n", cell.x, cell.y, cell.z, i, j, k);
                        result = false;
                        continue;
                    }
                    
                    cCellDelta3 scell = delta[cell.x + w * cell.y + w * h * cell.z];
                    
                    if (scell.x != 0 || scell.y != 0 || scell.z != 0)
                    {
                        result = false;
                        fprintf(s, "non-boundary cp (%d:%d:%d) @ cell %d:%d:%d\n", scell.x, scell.y, scell.z, i, j, k);
                    }
                }

        return result;
    }

    // Distance field testing

    void GenerateDistance(int kX, int kY, uint32_t mask[], int method, FILE* s)
    {
        cCellDelta2 delta[kX * kY];
        
        InitDeltaFromBitMask(kX, kY, mask, delta);

        if (s)
        {
            PrintDeltas(kX, kY, delta, s);
            fprintf(s, "\n");
        }

        switch (method)
        {
        case 0:
            if (s)
                fprintf(s, "\nDanielsson:\n");
            Danielsson(kX, kY, delta);
            break;
        case 1:
            if (s)
                fprintf(s, "\nFast Sweep:\n");
            FastSweep(kX, kY, delta);
            break;
        case 2:
            if (s)
                fprintf(s, "\nDanielsson x 2:\n");
            Danielsson(kX, kY, delta);
            Danielsson(kX, kY, delta);
            break;
        case 3:
            if (s)
                fprintf(s, "\nFast Sweep x 2:\n");
            FastSweep(kX, kY, delta);
            FastSweep(kX, kY, delta);
            break;
        default:
            break;
        }
        
        if (s)
        {
            PrintDeltas(kX, kY, delta, s);
            if (!CheckDeltas(kX, kY, delta, s))
                return;

            fprintf(s, "\nOptimal:\n");
            
            cCellDelta2 optDelta[kX * kY];
            InitDeltaFromBitMask(kX, kY, mask, optDelta);
            FindOptimalDeltas(kX, kY, optDelta);
            PrintDeltas(kX, kY, optDelta, s);
            
            fprintf(s, "\nDiff from optimal:\n");
            PrintDiff(kX, kY, delta, optDelta, s);
        }

    #if USE_TGA
        uint8_t dist8[kX * kY];
        float maxLen = sqrtf(kX * kX + kY * kY) * 0.125f;
        
        for (int i = 0, n = kX * kY; i < n; i++)
        {
            float dist2 = delta[i].x * delta[i].x + delta[i].y * delta[i].y;
            float normDist = sqrtf(dist2) / maxLen;
            if (normDist > 1.0f)
                normDist = 1.0f;
                
            dist8[i] = lrintf(normDist * 255);
        }
        
        tga_image tga;
        init_tga_image(&tga, dist8, kX, kY, 8);
        tga.image_type = TGA_IMAGE_TYPE_MONO;
        tga.image_descriptor &= ~TGA_T_TO_B_BIT;
        tga_write("dist.tga", &tga);
    #endif
    }

    void TestImageDistance(int w, int h, int seed, int method, FILE* s)
    {
        const int kSX = (w + 31) / 32;

        uint32_t mask[kSX * h];
        assert(MaskSize(w, h) == sizeof(mask) / sizeof(uint32_t));

        memset(mask, 0, sizeof(mask));

        if (seed < 0 && h >= 16)
            CreateBitMaskFromBlock(w, h, -seed, mask);
        else
            CreateBitMask(w, h, seed, mask);
        
        if (s)
            PrintMask(w, h, mask, s);

        GenerateDistance(w, h, mask, method, s);
    }

#if USE_TGA
    void TestImageDistanceFromFile(const char* filename, int method, FILE* s)
    {
        tga_image image;
        tga_result result = tga_read(&image, filename);

        if (result != TGA_NOERR)
        {
            printf("error reading %s: %s\n", filename, tga_error(result));
            return;
        }

        int w = image.width;
        int h = image.height;
        int pb = (image.pixel_depth + 7) / 8;
        const int sw = (w + 31) / 32;

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

                    uint8_t b, g, r, a;
                    tga_unpack_pixel(image.image_data + (y * w + x) * pb, image.pixel_depth, &b, &g, &r, &a);

                    if (r < 128)
                        *maskCell |= 1 << j;
                }
            }

        GenerateDistance(w, h, mask, method, s);
    }
#endif

    void TestVolumeDistance(const int kX, const int kY, const int kZ, int seed, int method, FILE* s)
    {
        const int kSX = (kX + 31) / 32;
        
        if (s)
            fprintf(s, "Init seed=%d:\n", seed);
        
        uint32_t mask[kSX * kY * kZ];
        assert(MaskSize(kX, kY, kZ) == sizeof(mask) / sizeof(uint32_t));

        memset(mask, 0, sizeof(mask));
        
        if (seed < 0)
            CreateBitMaskFromBlock(kX, kY, kZ, -seed, mask);
        else
            CreateBitMask(kX, kY, kZ, seed, mask);

        if (s)
            PrintMask(kX, kY, kZ, mask, s);
        
        cCellDelta3 delta[kX * kY * kZ];
        
        InitDeltaFromBitMask(kX, kY, kZ, mask, delta);

        if (s)
            fprintf(s, "\nMethod %d:\n", method);

        switch (method)
        {
        case 0:
            if (s)
                fprintf(s, "\nDanielsson:\n");
            Danielsson(kX, kY, kZ, delta);
            break;
        case 1:
            if (s)
                fprintf(s, "\nFast Sweep:\n");
            FastSweep(kX, kY, kZ, delta);
            break;
        case 2:
            if (s)
                fprintf(s, "\nFast Sweep x 4:\n");
            Danielsson(kX, kY, kZ, delta);
            Danielsson(kX, kY, kZ, delta);
            Danielsson(kX, kY, kZ, delta);
            Danielsson(kX, kY, kZ, delta);
            break;
        case 3:
            if (s)
                fprintf(s, "\nFast Sweep x 4:\n");
            FastSweep(kX, kY, kZ, delta);
            FastSweep(kX, kY, kZ, delta);
            FastSweep(kX, kY, kZ, delta);
            FastSweep(kX, kY, kZ, delta);
            break;
        default:
            if (s)
                fprintf(s, "No such method\n");
            break;
        }

        if (s)
        {
            PrintDeltas(kX, kY, kZ, delta, s);

            if (!CheckDeltas(kX, kY, kZ, delta, s))
                return;
        
            fprintf(s, "\nOptimal:\n");
        
            cCellDelta3 optDelta[kX * kY * kZ];
            InitDeltaFromBitMask(kX, kY, kZ, mask, optDelta);
            FindOptimalDeltas(kX, kY, kZ, optDelta);

            PrintDeltas(kX, kY, kZ, optDelta, s);

            fprintf(s, "\nDiff from optimal:\n");
            PrintDiff(kX, kY, kZ, delta, optDelta, s);
        }

    #if USE_TGA
        // Write as images
        uint8_t dist8[kX * kY];
        float maxLen = 10.0f;//sqrtf(kX * kX + kY * kY);
        
        for (int slice = 0; slice < kZ; slice++)
        {
            cCellDelta3* spos = delta + kX * kY * slice;
            
            for (int i = 0, n = kX * kY; i < n; i++)
            {
                float sx = spos[i].x;
                float sy = spos[i].y;
                float sz = spos[i].z;

                float dist2 = sx * sx + sy * sy + sz * sz;
                
                float normDist = sqrtf(dist2) / maxLen;
                
                if (normDist > 1.0f)
                    normDist = 1.0f;
                
                dist8[i] = lrintf(normDist * 255);
            }
        
            tga_image tga;
            init_tga_image(&tga, dist8, kX, kY, 8);
            tga.image_type = TGA_IMAGE_TYPE_MONO;
            tga.image_descriptor &= ~TGA_T_TO_B_BIT;
            
            char sliceTGAName[128];
            sprintf(sliceTGAName, "dist-%03d.tga", slice);
            tga_write(sliceTGAName, &tga);        
        }
    #endif
    }

    // Occlusion testing

    void GenerateOcclusion(const int kX, const int kY, float* dirW[4], int method, FILE* s)
    {
        if (method == 0)
            OcclusionSweep(kX, kY, dirW[0]);

        if (s)
            for (int i = 0; i < 4; i++)
            {
                fprintf(s, "quadrant %d:\n", i);
                PrintDistances(kX, kY, dirW[i], s);
            }

    #if USE_TGA
        uint8_t ao[kX * kY];

        for (int i = 0, n = kX * kY; i < n; i++)
        {
            float aoF = 1.0f - 0.25f * (dirW[0][i] + dirW[1][i] + dirW[2][i] + dirW[3][i]);
            
            if (aoF < 0.0f)
                aoF = 0.0f;
                
            ao[i] = lrintf(aoF * 255);
        }

        tga_image tga;
        init_tga_image(&tga, ao, kX, kY, 8);
        tga.image_type = TGA_IMAGE_TYPE_MONO;
        tga.image_descriptor &= ~TGA_T_TO_B_BIT;
        tga_write("AO.tga", &tga);    

        uint8_t aoDir[kX * kY][4];
        
        for (int i = 0, n = kX * kY; i < n; i++)
        {
            float aoFW = 0.25f * (dirW[0][i] + dirW[1][i] + dirW[2][i] + dirW[3][i]);
            float aoFX = dirW[0][i] - dirW[1][i] + dirW[2][i] - dirW[3][i];
            float aoFY = dirW[0][i] + dirW[1][i] - dirW[2][i] - dirW[3][i];
            
            float invLen = 1.0f / sqrtf(aoFX * aoFX + aoFY * aoFY);
            
            aoFX = aoFW * invLen * aoFX * 0.5f + 0.5f;
            aoFY = aoFW * invLen * aoFY * 0.5f + 0.5f;

            aoDir[i][2] = lrintf(aoFX * 255);
            aoDir[i][1] = lrintf(aoFY * 255);
            aoDir[i][0] = 128; // lrintf(aoFW * 255);
            aoDir[i][3] = 0;
        }
        
        init_tga_image(&tga, aoDir[0], kX, kY, 32);
        tga.image_type = TGA_IMAGE_TYPE_BGR;
        tga.image_descriptor &= ~TGA_T_TO_B_BIT;
        tga_write("AODir.tga", &tga);    

        for (int d = 0; d < 4; d++)
        {
            for (int i = 0, n = kX * kY; i < n; i++)
            {
                uint8_t c = lrintf((1.0f - dirW[d][i]) * 255);
                
                aoDir[i][2] = c;
                aoDir[i][1] = c;
                aoDir[i][0] = c;
                aoDir[i][3] = 0xFF;
            }
        
            init_tga_image(&tga, aoDir[0], kX, kY, 32);
            tga.image_type = TGA_IMAGE_TYPE_BGR;
            tga.image_descriptor &= ~TGA_T_TO_B_BIT;
            char cmptName[128];
            sprintf(cmptName, "AOC-%d.tga", d);
            tga_write(cmptName, &tga);    
        }
    #endif
    }

    void TestImageOcclusion(const int kX, const int kY, int seed, int method, FILE* s)
    {
        const int kSX = (kX + 31) / 32;

        uint32_t mask[kSX * kY];
        assert(MaskSize(kX, kY) == sizeof(mask) / sizeof(uint32_t));

        memset(mask, 0, sizeof(mask));

        if (seed < 0)
            CreateBitMaskFromBlock(kX, kY, -seed, mask);
        else
            CreateBitMask(kX, kY, seed, mask);

        float dirW[4][kX * kY];
        InitDirWFromBitMask(kX, kY, mask, dirW[0]);

        float* dirW4[4] = { dirW[0], dirW[1], dirW[2], dirW[3] };
        GenerateOcclusion(kX, kY, dirW4, method, s);
    }

#if USE_TGA
    void TestImageOcclusionFromFile(const char* filename, int method, FILE* s)
    {
        tga_image image;
        tga_result result = tga_read(&image, filename);

        if (result != TGA_NOERR)
        {
            printf("error reading %s: %s\n", filename, tga_error(result));
            return;
        }

        int w = image.width;
        int h = image.height;
        int n = w * h;
        int pb = (image.pixel_depth + 7) / 8;

        float* dirW4[4];

        for (int i = 0; i < 4; i++)
            dirW4[i] = new float[n];

        for (int i = 0; i < n; i++)
        {
            uint8_t b, g, r, a;
            tga_unpack_pixel(image.image_data + i * pb, image.pixel_depth, &b, &g, &r, &a);

            dirW4[0][i] = 1.0f - r * (1.0f / 255.0f);
        }

        for (int i = 1; i < 4; i++)
            memcpy(dirW4[i], dirW4[0], sizeof(float) * w * h);

        GenerateOcclusion(w, h, dirW4, method, s);
    }
#endif

    void GenerateOcclusion(const int kX, const int kY, const int kZ, float* dirW, int method, FILE* s)
    {
        const size_t volumeStride = kX * kY * kZ;

        if (s)
            fprintf(s, "\nMethod %d:\n", method);

        if (method == 0)
            OcclusionSweep(kX, kY, kZ, dirW);
        else if (method == 1)
        {
            size_t volumeOctantsSize = sizeof(float) * 8 * volumeStride;
            float* outDirW = (float*) malloc(volumeOctantsSize);

            OcclusionSweep(kX, kY, kZ, dirW, outDirW);
            memcpy(dirW, outDirW, volumeOctantsSize);

            free(outDirW);
        }

        if (s)
            for (int i = 0; i < 8; i++)
            {
                fprintf(s, "octant %d:\n", i);
                PrintDistances(kX, kY, kZ, dirW + volumeStride * i, s);
            }

    #if USE_TGA
        // Write as images
        const size_t sliceStride = kX * kY;
        uint8_t ao[kX * kY];
        uint8_t aoDir[kX * kY][4];
        
        for (int slice = 0; slice < kZ; slice++)
        {
            float* sDirW[8];

            for (int d = 0; d < 8; d++)
                sDirW[d] = dirW + sliceStride * slice + volumeStride * d;

            for (int i = 0, n = kX * kY; i < n; i++)
            {
                float aoF = sDirW[0][i];
                
                for (int d = 1; d < 8; d++)
                    aoF += sDirW[d][i];

                aoF = 1.0f - 0.125f * aoF;
                if (aoF < 0.0f)
                    aoF = 0.0f;

                ao[i] = lrintf(aoF * 255);
            }
        
            tga_image tga;
            init_tga_image(&tga, ao, kX, kY, 8);
            tga.image_type = TGA_IMAGE_TYPE_MONO;
            tga.image_descriptor &= ~TGA_T_TO_B_BIT;
            
            char sliceTGAName[128];
            sprintf(sliceTGAName, "AOV-%03d.tga", slice);
            tga_write(sliceTGAName, &tga);
            
            for (int i = 0, n = kX * kY; i < n; i++)
            {
                float aoFW = sDirW[0][i] + sDirW[1][i] + sDirW[2][i] + sDirW[3][i]
                           + sDirW[4][i] + sDirW[5][i] + sDirW[6][i] + sDirW[7][i];

                aoFW *= 0.125f;
                                    
                float aoFX = sDirW[0][i] - sDirW[1][i] + sDirW[2][i] - sDirW[3][i]
                           + sDirW[4][i] - sDirW[5][i] + sDirW[6][i] - sDirW[7][i];
                float aoFY = sDirW[0][i] + sDirW[1][i] - sDirW[2][i] - sDirW[3][i]
                           + sDirW[4][i] + sDirW[5][i] - sDirW[6][i] - sDirW[7][i];
                float aoFZ = sDirW[0][i] + sDirW[1][i] + sDirW[2][i] + sDirW[3][i]
                           - sDirW[4][i] - sDirW[5][i] - sDirW[6][i] - sDirW[7][i];
                
                float invLen = 1.0f / sqrtf(aoFX * aoFX + aoFY * aoFY + aoFZ * aoFZ);
                
                aoFX = aoFW * invLen * aoFX * 0.5f + 0.5f;
                aoFY = aoFW * invLen * aoFY * 0.5f + 0.5f;
                aoFZ = aoFW * invLen * aoFZ * 0.5f + 0.5f;
                
                aoDir[i][2] = lrintf(aoFX * 255);
                aoDir[i][1] = lrintf(aoFY * 255);
                aoDir[i][0] = lrintf(aoFZ * 255);
                aoDir[i][3] = 0; // lrintf(aoFW * 255);
            }
            
            init_tga_image(&tga, aoDir[0], kX, kY, 32);
            tga.image_type = TGA_IMAGE_TYPE_BGR;
            tga.image_descriptor &= ~TGA_T_TO_B_BIT;
            
            sprintf(sliceTGAName, "AOVDir-%03d.tga", slice);
            tga_write(sliceTGAName, &tga);
        }
    #endif
    }
    
    void TestVolumeOcclusion(const int kX, const int kY, const int kZ, int seed, int method, FILE* s)
    {
        const size_t volumeStride = kX * kY * kZ;
        const size_t kSX = (kX + 31) / 32;

        uint32_t mask[kSX * kY * kZ];
        assert(MaskSize(kX, kY, kZ) == sizeof(mask) / sizeof(uint32_t));

        memset(mask, 0, sizeof(mask));

        if (seed < 0)
            CreateBitMaskFromBlock(kX, kY, kZ, -seed, mask);
        else
            CreateBitMask(kX, kY, kZ, seed, mask);

        // Switch to malloc as volume data can be large.
        size_t volumeOctantsSize = sizeof(float) * 8 * volumeStride;
        float* dirW = (float*) malloc(volumeOctantsSize);
        InitDirWFromBitMask(kX, kY, kZ, mask, dirW);

        GenerateOcclusion(kX, kY, kZ, dirW, method, s);

        free(dirW);
    }

#if USE_OBJ
    void TestVolumeOcclusionFromFile(const int w, const int h, const int d, const char* sourceFile, int method, FILE* s)
    {
        FILE* meshFile = fopen(sourceFile, "r");

        if (!meshFile)
        {
            fprintf(stderr, "Can't open %s\n", sourceFile);
            return;
        }

        cMesh mesh;
        ReadObjFile(meshFile, &mesh);

        const size_t volumeStride = w * h * d;

        // Switch to malloc as volume data can be large.
        size_t volumeOctantsSize = sizeof(float) * 8 * volumeStride;
        float* dirW = (float*) malloc(volumeOctantsSize);

        Bounds3f bbox = { { +FLT_MAX, +FLT_MAX, +FLT_MAX }, { -FLT_MAX, -FLT_MAX, -FLT_MAX } };

        for (int i = 0, n = mesh.mPositions.size(); i < n; i++)
        {
            Vec3f p = mesh.mPositions[i];

            if (bbox.mMin.x > p.x)
                bbox.mMin.x = p.x;
            else
            if (bbox.mMax.x < p.x)
                bbox.mMax.x = p.x;
            
            if (bbox.mMin.y > p.y)
                bbox.mMin.y = p.y;
            else
            if (bbox.mMax.y < p.y)
                bbox.mMax.y = p.y;
            
            if (bbox.mMin.z > p.z)
                bbox.mMin.z = p.z;
            else
            if (bbox.mMax.z < p.z)
                bbox.mMax.z = p.z;
        }

        bbox = FindAOBounds(0.5f, bbox);

        bool result = CreateDirW8FromTriangles
        (
            mesh.mPositionIndices.size() / 3,
            mesh.mPositionIndices.data(),
            mesh.mPositions.data(),
            bbox,
            w, h, d,
            dirW
        );

        if (result)
            GenerateOcclusion(w, h, d, dirW, method, s);

        free(dirW);
    }
#endif
}

int main(int argc, const char* argv[])
{
    if (argc == 1)
    {
        printf
        (
            "Distance/Occlusion experimental tool\n"
            "\n"
            "Options:\n"
            "  -d: generate image distance field. Methods 0: Danielsson, 1: Sweep, 2: 2 x Danielsson, 3: 2 x Sweep\n"
            "  -D: generate volume distance field. Methods 0: Danielsson, 1: Sweep, 2: 4 x Danielsson, 3: 4 x Sweep\n"
            "  -o: generate image occlusion field\n"
            "  -O: generate volume occlusion field. Methods 0: standard, 1: no self-occlusion\n"
            "  -m <int>: select method variant\n"
            "\n"
            "  -s w [h [d]]: set dimensions of source\n"
            "  -x <int>: use random point image/volume with the given seed\n"
            "  -b <int>: use test block image/volume with the given number of sides\n"
        #if USE_TGA
            "  -f <filename>: use tga as input image\n"
        #endif
        #if USE_OBJ
            "  -f <filename>: use obj to define input volume\n"
        #endif
        #if USE_TGA || USE_OBJ
            "\n"
            "  -l : log detailed output\n"
        #endif
        );
        
        return -1;
    }

    argc--;
    argv++;

    enum tAlgorithm
    {
        kImageDistance,
        kVolumeDistance,
        kImageOcclusion,
        kVolumeOcclusion,
        kMaxAlgorithms
    };

    int algorithm = 0;
    int method = 0;
    int dim[3] = { 0, 0, 0 };
    int seed = 1;
#if USE_TGA || USE_OBJ
    const char* sourceFile = 0;
    FILE* out = 0;
#else
    FILE* out = stdout;
#endif

    while (argc > 0 && argv[0][0] == '-')
    {
        const char* option = argv[0] + 1;

        argc--;
        argv++;

        switch (option[0])
        {
        case 'd':
            algorithm = kImageDistance;
            break;
        case 'D':
            algorithm = kVolumeDistance;
            break;
        case 'o':
            algorithm = kImageOcclusion;
            break;
        case 'O':
            algorithm = kVolumeOcclusion;
            break;
        case 'l':
            out = stdout;
            break;
        case 'm':
            if (argc > 0)
            {
                method = atoi(argv[0]);

                argc--;
                argv++;
            }
            break;
        case 's':
            for (int i = 0; i < 3 && argc > 0 && argv[0][0] != '-'; i++)
            {
                dim[i] = atoi(argv[0]);

                argc--;
                argv++;
            }
            break;
        case 'x':
            if (argc > 0)
            {
                seed = atoi(argv[0]);

                argc--;
                argv++;
            }
            break;
        case 'b':
            if (argc > 0)
            {
                seed = -atoi(argv[0]);

                argc--;
                argv++;
            }
            break;
    #if USE_TGA || USE_OBJ
        case 'f':
            if (argc > 0)
            {
                sourceFile = argv[0];

                argc--;
                argv++;
            }
            break;
    #endif
        default:
            printf("unknown option: %s\n", option);
        }
    }

    if (dim[0] == 0)
        dim[0] = (algorithm == kVolumeDistance) || (algorithm == kVolumeOcclusion) ? 16 : 64;

    if (dim[1] == 0)
        dim[1] = dim[0];

    if (dim[2] == 0)
        dim[2] = dim[1];

    switch (algorithm)
    {
    case kImageDistance:
    #if USE_TGA
        if (sourceFile)
            TestImageDistanceFromFile(sourceFile, method, out);
        else
    #endif
        TestImageDistance(dim[0], dim[1], seed, method, out);
        break;

    case kVolumeDistance:
    #if USE_OBJ
// XXX        if (sourceFile)
//            TestVolumeDistanceFromFile(dim[0], dim[1], dim[2], sourceFile, method, out);
//        else
    #endif
        TestVolumeDistance(dim[0], dim[1], dim[2], seed, method, out);
        break;

    case kImageOcclusion:
    #if USE_TGA
        if (sourceFile)
            TestImageOcclusionFromFile(sourceFile, method, out);
        else
    #endif
        TestImageOcclusion(dim[0], dim[1], seed, method, out);
        break;

    case kVolumeOcclusion:
    #if USE_OBJ
        if (sourceFile)
            TestVolumeOcclusionFromFile(dim[0], dim[1], dim[2], sourceFile, method, out);
        else
    #endif
        TestVolumeOcclusion(dim[0], dim[1], dim[2], seed, method, out);
        break;

    default:
        printf("unknown algorithm\n");
    }
}
