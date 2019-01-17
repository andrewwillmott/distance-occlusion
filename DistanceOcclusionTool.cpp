//------------------------------------------------------------------------------
// Purpose: Distance Field and Occlusion generator tool
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#include "DistanceOcclusionLib.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
    #pragma warning (disable: 4996)     // let's not make fopen an error.
    #define strlcpy(d, s, ds) strcpy_s(d, ds, s)
#elif defined(__GNUC__) && !defined(__clang__)
    #define strlcpy strncpy // thanks glibc guys
#endif

#ifndef USE_STB
    #define USE_STB 1
#endif

#ifndef USE_MESH
    #define USE_MESH 1
#endif

#if USE_STB
    #include "stb_image_mini.h"
#endif

#if USE_MESH
    #include "MeshSupport.cpp"
    using namespace MSL;
#endif

using namespace DOL;

namespace
{
    ////////////////////////////////////////////////////////////////////////////
    // MARK: - Utilities
    ////////////////////////////////////////////////////////////////////////////

    inline int32_t Dist2(const cCellDelta2& cd)
    {
        return cd.x * cd.x + cd.y * cd.y;
    }
    inline int32_t Dist2(const cCellDelta3& cd)
    {
        return cd.x * cd.x + cd.y * cd.y + cd.z * cd.z;
    }
    inline float Dist(const cCellDelta2& cd)
    {
        return sqrtf((float) Dist2(cd));
    }
    inline float Dist(const cCellDelta3& cd)
    {
        return sqrtf((float) Dist2(cd));
    }
    inline bool IsBoundary(const cCellDelta2& cd, int cw = 1)
    {
        return abs(cd.x) < cw && abs(cd.y) < cw;
    }
    inline bool IsBoundary(const cCellDelta3& cd, int cw = 1)
    {
        return abs(cd.x) < cw && abs(cd.y) < cw && abs(cd.z) < cw;
    }

    inline uint8_t EncodeU8Signed(float v)
    {
        v = (v < -1.0f) ? -1.0f : ((v > 1.0f) ? 1.0f : v);
        return (uint8_t) lrintf(128 + v * 127); // ensure 0 = 128
    }

    // MARK: - Printing

    void PrintMask(int w, int h, const uint32_t* mask, FILE* s)
    {
        for (int y = 0; y < h; y++)
        {
            for (int j = 0; j < w; j += 32)
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
        for (int z = 0; z < d; z++)
        {
            for (int y = 0; y < h; y++)
            {
                for (int j = 0; j < w; j += 32)
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

    void PrintDistances(int w, int h, const int16_t distances[], FILE* s)
    {
        for (int i = 0; i < h; i++)
        {
            for (int i = 0; i < w - 1; i++)
                fprintf(s, "%2d ", (*distances++));
            
            fprintf(s, "%2d\n", (*distances++));
        }
    }

    void PrintDistances(int w, int h, const float* dist, FILE* s)
    {    
        for (int i = 0; i < h; i++)
        {
            for (int i = 0; i < w - 1; i++)
                fprintf(s, "%1.3f ", (*dist++));
            
            fprintf(s, "%1.3f\n", (*dist++));
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
                    fprintf(s, "%1.3f ", (*dist++));
                    
                fprintf(s, "%1.3f\n", (*dist++));
            }
        }
    }

    void PrintDeltas(int w, int h, const cCellDelta2 deltas[], FILE* s)
    {    
        for (int j = 0; j < h; j++)
            for (int i = 0; i < w; i++, deltas++)
            {
                if (deltas->x >= kMaxDelta2)
                    fprintf(s, "   x   ");
                else
                    fprintf(s, "%3d %3d", deltas->x, deltas->y);

                if (i < w - 1)
                    fprintf(s, "|");
                else
                    fprintf(s, "\n");
            }

        fprintf(s, "\n");
    }

    void PrintDeltas(int w, int h, int d, const cCellDelta3 deltas[], FILE* s)
    {    
        for (int k = 0; k < d; k++)
        {
            fprintf(s, "slice %d>\n", k);
            
            for (int j = 0; j < h; j++)
                for (int i = 0; i < w; i++, deltas++)
                {
                    if (deltas->x >= kMaxDelta3)
                        fprintf(s, "     x     ");
                    else
                        fprintf(s, "%3d %3d %3d", deltas->x, deltas->y, deltas->z);

                    if (i < w - 1)
                        fprintf(s, "|");
                    else
                        fprintf(s, "\n");
                }
        }
    }

    void PrintDiff(int w, int h, const cCellDelta2 deltas[], const cCellDelta2 refs[], FILE* s)
    {
        for (int i = 0; i < h; i++)
        {
            for (int j = 0; j < w; j++, deltas++, refs++)
            {
                if (deltas->x == refs->x && deltas->y == refs->y)
                    fprintf(s, "   -   ");
                else
                {
                    float dist = Dist(*deltas) - Dist(*refs);
                    
                    if (fabsf(dist) < 1e-6f)
                        fprintf(s, "   x   ");
                    else
                        fprintf(s, "%7.2f", dist);
                }
                
                if (j < w - 1)
                    fprintf(s, ", ");
                else
                    fprintf(s, "\n");
            }
        }
    }

    void PrintDiff(int w, int h, const float distances[], const cCellDelta2 refs[], FILE* s)
    {
        for (int i = 0; i < h; i++)
        {
            for (int j = 0; j < w; j++, distances++, refs++)
            {
                float dist = *distances - Dist(*refs);
                if (dist == 0.0f)
                    fprintf(s, "   -   ");
                else if (fabsf(dist) < 1e-6f)
                        fprintf(s, "   x   ");
                else
                    fprintf(s, "%7.2f", dist);
                
                if (j < w - 1)
                    fprintf(s, ", ");
                else
                    fprintf(s, "\n");
            }
        }
    }

    bool PrintError(int w, int h, const cCellDelta2 deltas[], const cCellDelta2 refs[], FILE* s)
    {
        float maxDist = 0.0f;

        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j++, deltas++, refs++)
            {
                float dist = Dist(*deltas) - Dist(*refs);

                if (maxDist < dist)
                    maxDist = dist;
            }
            
        fprintf(s, "Max error: %g\n", maxDist);

        if (maxDist >= 1.0f)
        {
            fprintf(s, "FAILED\n");
            return false;
        }
        
        return true;
    }

    bool PrintError(int w, int h, const float distances[], const cCellDelta2 refs[], FILE* s)
    {
        float maxDist = 0.0f;

        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j++, distances++, refs++)
            {
                float dist = *distances - Dist(*refs);

                if (maxDist < dist)
                    maxDist = dist;
            }
        
        fprintf(s, "Max error: %g\n", maxDist);

        if (maxDist >= 1.0f)
        {
            fprintf(s, "FAILED\n");
            return false;
        }
        
        return true;
    }

    void PrintDiff(int w, int h, int d, const cCellDelta3 deltas[], const cCellDelta3* refs, FILE* s)
    {
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
                        float dist = Dist(*deltas) - Dist(*refs);

                        if (fabsf(dist) < 1e-6f)
                            fprintf(s, "   x   ");
                        else
                            fprintf(s, "%7.2f", dist);
                    }
                    
                    if (j < w - 1)
                        fprintf(s, ", ");
                    else
                        fprintf(s, "\n");
                }
            }
        }
    }

    bool PrintError(int w, int h, int d, const cCellDelta3 deltas[], const cCellDelta3* refs, FILE* s)
    {
        float maxDist = 0.0f;

        for (int k = 0; k < d; k++)
        for (int i = 0; i < h; i++)
            for (int j = 0; j < w; j++, deltas++, refs++)
            {
                float dist = Dist(*deltas) - Dist(*refs);

                if (maxDist < dist)
                    maxDist = dist;
            }
        
        fprintf(s, "Max error: %g\n", maxDist);

        if (maxDist >= 1.0f)
        {
            fprintf(s, "FAILED\n");
            return false;
        }
        
        return true;
    }

    // MARK: - Sanity checking

    bool CheckDeltas(int w, int h, cCellDelta2* deltas, FILE* s, int cw = 1)
    {
        bool result = true;
        
        for (int j = 0; j < h; j++)
            for (int i = 0; i < w; i++)
            {
                cCellDelta2 cell = deltas[i + w * j];
                cell.x += i * cw;
                cell.y += j * cw;
                
                if (cell.x < 0 || cell.x >= cw * w || cell.y < 0 || cell.y >= cw * h)
                {
                    fprintf(s, "out of bounds delta (%d:%d) @ cell %d:%d\n", cell.x, cell.y, i, j);
                    result = false;
                    continue;
                }
                
                int ri = cell.x / cw;
                int rj = cell.y / cw;

                cCellDelta2 rcell = deltas[ri + w * rj];
                
                if (!IsBoundary(rcell, cw))
                {
                    result = false;
                    fprintf(s, "non-boundary delta (%d:%d) @ cell %d:%d referenced by %d:%d\n", rcell.x, rcell.y, ri, rj, i, j);
                }
            }
            
        return result;
    }

    bool CheckDeltas(int w, int h, int d, cCellDelta3* deltas, FILE* s, int cw = 1)
    {
        bool result = true;
        
        for (int k = 0; k < d; k++)
            for (int j = 0; j < h; j++)
                for (int i = 0; i < w; i++)
                {
                    cCellDelta3 cell = deltas[i + w * j + w * h * k];
                    
                    cell.x += i * cw;
                    cell.y += j * cw;
                    cell.z += k * cw;
                    
                    if (cell.x < 0 || cell.x >= cw * w || cell.y < 0 || cell.y >= cw * h || cell.z < 0 || cell.z >= cw * d)
                    {
                        fprintf(s, "out of bounds delta (%d:%d:%d) @ cell %d:%d:%d\n", cell.x, cell.y, cell.z, i, j, k);
                        result = false;
                        continue;
                    }
                    
                    int ri = cell.x / cw;
                    int rj = cell.y / cw;
                    int rk = cell.z / cw;

                    cCellDelta3 rcell = deltas[ri + w * rj + w * h * rk];

                    if (!IsBoundary(rcell, cw))
                    {
                        result = false;
                        fprintf(s, "non-boundary delta (%d:%d:%d) @ cell %d:%d:%d referenced by %d:%d:%d\n", rcell.x, rcell.y, rcell.z, ri, rj, rk, i, j, k);
                    }
                }

        return result;
    }


    // Check border-style init
    inline bool Occupied(const uint32_t mask[], int w, int h, int x, int y)
    {
        int ms = MaskSize(w);

        if (x < 0)
            x = 0;
        else if (x >= w)
            x = w - 1;

        if (y < 0)
            y = 0;
        else if (y >= h)
            y = h - 1;

        return (mask[y * ms + (x >> 5)] & (1 << (x & 0x1F))) != 0;
    }

    inline bool Expect(bool value, const uint32_t mask[], int w, int h, int x, int y, int dx, int dy, FILE* s)
    {
        if (Occupied(mask, w, h, x + dx, y + dy) != value)
        {
            fprintf(s, "bad cell (%d %d), offset (%d, %d) != %d\n", x, y, dx, dy, value);
            return false;
        }
        return true;
    }

    bool CheckBorderDeltas(int w, int h, const uint32_t mask[], cCellDelta2 deltas[], FILE* s)
    {
        bool good = true;

        for (int x = 0; x < h; x++)
            for (int y = 0; y < w; y++)
            {
                cCellDelta2* cell = deltas + y * w + x;
                bool occ = Occupied(mask, w, h, x, y);

                if (IsBoundary(*cell, 2))
                {
                    good = Expect(!occ, mask, w, h, x, y, cell->x, cell->y, s) && good;

                    if (cell->x != 0 && cell->y != 0)   // if diagonal check if there was a closer edge...
                    {
                        good = Expect(occ, mask, w, h, x, y, -1, 0, s) && good;
                        good = Expect(occ, mask, w, h, x, y, +1, 0, s) && good;
                        good = Expect(occ, mask, w, h, x, y, 0, -1, s) && good;
                        good = Expect(occ, mask, w, h, x, y, 0, +1, s) && good;
                    }
                }
                else
                {
                    for (int dy = -1; dy <= +1; dy++)
                    for (int dx = -1; dx <= +1; dx++)
                        good = Expect(occ, mask, w, h, x, y, dx, dy, s) && good;
                }
            }

        return good;
    }

    inline bool Occupied(const uint32_t mask[], int w, int h, int , int x, int y, int z)
    {
        int ms = MaskSize(w);

        return (mask[(z * h + y) * ms + (x >> 5)] & (1 << (x & 0x1F))) != 0;
    }

    inline bool Expect(const uint32_t mask[], int w, int h, int d, int x, int y, int z, bool value)
    {
        if (x < 0 || x >= w)
            return true;

        if (y < 0 || y >= h)
            return true;

        if (z < 0 || z >= d)
            return true;

        return Occupied(mask, w, h, d, x, y, z) == value;
    }

    inline bool Expect(bool value, const uint32_t mask[], int w, int h, int d, int x, int y, int z, int dx, int dy, int dz, FILE* s)
    {
        if (!Expect(mask, w, h, d, x + dx, y + dy, z + dz, value))
        {
            if (s)
                fprintf(s, "bad cell (%d %d %d), offset (%d, %d, %d) != %d\n", x, y, z, dx, dy, dz, value);
            return false;
        }
        return true;
    }

    bool CheckBorderDeltas(int w, int h, int d, const uint32_t mask[], cCellDelta3 deltas[], FILE* s)
    {
        bool good = true;

        for (int z = 0; z < d; z++)
        for (int x = 0; x < h; x++)
        for (int y = 0; y < w; y++)
        {
            cCellDelta3* cell = deltas + z * w * h + y * w + x;
            bool occ = Occupied(mask, w, h, d, x, y, z);

            if (IsBoundary(*cell, 2))
            {
                good = Expect(!occ, mask, w, h, d, x, y, z, cell->x, cell->y, cell->z, s) && good;

                if (cell->x != 0 && cell->y != 0 && cell->z != 0)   // if box diagonal check if there was a closer face diagonal...
                {
                    good = Expect(occ, mask, w, h, d, x, y, z, -1, -1, 0, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, -1, +1, 0, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, +1, -1, 0, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, +1, +1, 0, s) && good;

                    good = Expect(occ, mask, w, h, d, x, y, z, -1, 0, -1, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, -1, 0, +1, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, +1, 0, -1, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, +1, 0, +1, s) && good;

                    good = Expect(occ, mask, w, h, d, x, y, z, 0, -1, -1, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, 0, -1, +1, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, 0, +1, -1, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, 0, +1, +1, s) && good;
                }
                if (cell->x != 0 && cell->y != 0)   // if diagonal check if there was a closer edge...
                {
                    good = Expect(occ, mask, w, h, d, x, y, z, -1, 0, 0, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, +1, 0, 0, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, 0, -1, 0, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, 0, +1, 0, s) && good;
                }
                if (cell->y != 0 && cell->z != 0)   // if diagonal check if there was a closer edge...
                {
                    good = Expect(occ, mask, w, h, d, x, y, z, 0, -1, 0, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, 0, +1, 0, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, 0, 0, -1, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, 0, 0, +1, s) && good;
                }
                if (cell->z != 0 && cell->x != 0)   // if diagonal check if there was a closer edge...
                {
                    good = Expect(occ, mask, w, h, d, x, y, z, -1, 0,  0, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z, +1, 0,  0, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z,  0, 0, -1, s) && good;
                    good = Expect(occ, mask, w, h, d, x, y, z,  0, 0, +1, s) && good;
                }
            }
            else
            {
                for (int dz = -1; dz <= +1; dz++)
                for (int dy = -1; dy <= +1; dy++)
                for (int dx = -1; dx <= +1; dx++)
                    good = Expect(occ, mask, w, h, d, x, y, z, dx, dy, dz, s) && good;
            }
        }

        return good;
    }


    // MARK: - Image IO
    void WriteImage(const char* name, const char* type, int w, int h, uint8_t data[], int index = -1)
    {
    #if USE_STB
        char nameAndExt[128];

        if (index >= 0)
            sprintf(nameAndExt, "%s-%s-%03d.png", name, type, index);
        else
            sprintf(nameAndExt, "%s-%s.png", name, type);

        stbi_write_png(nameAndExt, w, h, 1, data, 0);
        printf("written %s\n", nameAndExt);
    #endif
    }

    void WriteImage(const char* name, const char* type, int w, int h, uint8_t data[][4], int index = -1)
    {
    #if USE_STB
        char nameAndExt[128];

        if (index >= 0)
            sprintf(nameAndExt, "%s-%s-%03d.png", name, type, index);
        else
            sprintf(nameAndExt, "%s-%s.png", name, type);

        stbi_write_png(nameAndExt, w, h, 4, data, 0);
        printf("written %s\n", nameAndExt);
    #endif
    }

#if USE_STB
    uint32_t* BitMaskFromImageFile(const char* filename, int threshold, bool lessThan, int& w, int& h)
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

                    uint8_t v;
                    v = data[y * w + x];

                    if ((v >= threshold) ^ lessThan)
                        *maskCell |= 1 << j;
                }
        }

        stbi_image_free(data);

        return mask;
    }
#endif


    ////////////////////////////////////////////////////////////////////////////
    // MARK: - Distance generation
    ////////////////////////////////////////////////////////////////////////////

    void WriteDistances(const char* name, const char* type, int w, int h, float maxLen, const cCellDelta2 deltas[])
    {
        uint8_t* dist8 = new uint8_t[w * h];

        for (int i = 0, n = w * h; i < n; i++)
            dist8[i] = EncodeU8(Dist(deltas[i]) / maxLen);

        WriteImage(name, type, w, h, dist8);
        delete[] dist8;
    }

    template<class T> void WriteDistances(const char* name, const char* type, int w, int h, float maxLen, const T distances[])
    {
        uint8_t* dist8 = new uint8_t[w * h];

        for (int i = 0, n = w * h; i < n; i++)
            dist8[i] = EncodeU8(distances[i] / maxLen);

        WriteImage(name, type, w, h, dist8);
        delete[] dist8;
    }

    enum tMethod
    {
        kMethodDanielsson,      // O(n), fastest delta technique, decent accuracy
        kMethodFastSweep,       // O(n), ~50% slower, close to perfect accuracy
        kMethodJumpFlood,       // O(n log n), GPU-oriented algorithm, for checking accuracy only
        kMethodBruteForce,      // O(n^2), accurate, for checking other algorithms at low n.
        kMethodChamfer,         // O(n), has issues on diagonals. Cheapest but quality suffers.
        kMethodFelzenszwalb,    // O(n) but much larger overhead than the others. Accurate.
        kNumMethods
    };
    
    void FindDeltas(int method, int w, int h, cCellDelta2 deltas[], FILE* s, int cw = 1)
    {
        switch (method)
        {
        case kMethodDanielsson:
            if (s)
                fprintf(s, "\nDanielsson:\n");
            Danielsson(w, h, deltas, cw);
            break;
        case kMethodFastSweep:
            if (s)
                fprintf(s, "\nFast Sweep:\n");
            FastSweep(w, h, deltas, cw);
            break;
        case kMethodJumpFlood:
            if (s)
                fprintf(s, "\nJump Flood:\n");
            JumpFlood(w, h, deltas, cw);
            break;
        case kMethodBruteForce:
            if (s)
                fprintf(s, "\nBrute Force:\n");
            BruteForce(w, h, deltas, cw);
            break;
        default:
            fprintf(stderr, "\nUnknown method for this technique\n");
            break;
        }
    }

    bool GenerateDistanceFelzenszwalb(const char* name, int w, int h, float maxLen, uint32_t mask[], FILE* s, bool check)
    {
        // This does not use the more accurate (and more generalisable) delta cell approach, so has a different flow
        float* distances = new float[w * h];    // one of these things...

        InitDistancesFromBitMask(w, h, mask, distances);

        if (s)
            PrintDistances(w, h, distances, s);

        Felzenszwalb(w, h, distances);

        for (int i = 0; i < w * h; i++)
            distances[i] = sqrtf(distances[i]);

        if (s)
            PrintDistances(w, h, distances, s);
        
        if (check)
        {
            cCellDelta2* optDeltas = new cCellDelta2[w * h];
            InitDeltasFromBitMask(w, h, mask, optDeltas);
            BruteForce(w, h, optDeltas);

            if (s)
            {
                fprintf(s, "\nOptimal:\n");
                PrintDeltas(w, h, optDeltas, s);
                fprintf(s, "Diff from optimal:\n");
                PrintDiff(w, h, distances, optDeltas, s);
            }
            
            PrintError(w, h, distances, optDeltas, stderr);
            
            delete[] optDeltas;
        }

        WriteDistances(name, "df", w, h, maxLen, distances);

        delete[] distances;
        return true;
    }

    bool GenerateDistanceChamfer(const char* name, int w, int h, float maxLen, uint32_t mask[], FILE* s, bool check)
    {
        // This does not use the more accurate (and more generalisable) delta cell approach, so has a different flow
        int16_t* distances = new int16_t[w * h];    // one of these things...

        InitDistancesFromBitMask(w, h, mask, distances, INT16_MAX - 16);    // -16 to avoid overflow with Chamfer adds

        if (s)
            PrintDistances(w, h, distances, s);

        int cm = Chamfer(w, h, distances);

        if (s)
            PrintDistances(w, h, distances, s);
        
        if (check)
        {
            cCellDelta2* deltas = new cCellDelta2[w * h];   // easiest just to convert to delta form
            for (int i = 0; i < w * h; i++)
            {
                deltas[i].x = distances[i];
                deltas[i].y = 0;
            }

            cCellDelta2* optDeltas = new cCellDelta2[w * h];
            InitDeltasFromBitMask(w, h, mask, optDeltas);
            BruteForce(w, h, optDeltas, cm);

            if (s)
            {
                fprintf(s, "\nOptimal:\n");
                PrintDeltas(w, h, optDeltas, s);
                fprintf(s, "Diff from optimal:\n");
                PrintDiff(w, h, deltas, optDeltas, s);
            }
            
            PrintError(w, h, deltas, optDeltas, stderr);
            
            delete[] deltas;
            delete[] optDeltas;
        }

        WriteDistances(name, "df", w, h, maxLen * cm, distances);
        
        delete[] distances;
        return true;
    }

    bool GenerateDistance(const char* name, int w, int h, float maxLen, uint32_t mask[], int method, FILE* s, bool check)
    {
        // the non-delta approaches have special cases
        if (method == kMethodChamfer)
            return GenerateDistanceChamfer(name, w, h, maxLen, mask, s, check);
        else if (method == kMethodFelzenszwalb)
            return GenerateDistanceFelzenszwalb(name, w, h, maxLen, mask, s, check);
    
        cCellDelta2* deltas = new cCellDelta2[w * h];
        
        InitDeltasFromBitMask(w, h, mask, deltas);

        if (s)
        {
            PrintDeltas(w, h, deltas, s);
            fprintf(s, "Generating exterior distances\n");
        }

        FindDeltas(method, w, h, deltas, s);

        if (s)
            PrintDeltas(w, h, deltas, s);

        if (check)
        {
            if (!CheckDeltas(w, h, deltas, stderr))
                return false;

            if (method != kMethodBruteForce)
            {
                cCellDelta2* optDeltas = new cCellDelta2[w * h];
                InitDeltasFromBitMask(w, h, mask, optDeltas);
                BruteForce(w, h, optDeltas);

                if (s)
                {
                    fprintf(s, "\nOptimal:\n");
                    PrintDeltas(w, h, optDeltas, s);
                    fprintf(s, "Diff from optimal:\n");
                    PrintDiff(w, h, deltas, optDeltas, s);
                }
                
                PrintError(w, h, deltas, optDeltas, stderr);
                delete[] optDeltas;
            }
        }

        WriteDistances(name, "df", w, h, maxLen, deltas);
        delete[] deltas;

        return true;
    }

    // Volumes
    void WriteDistances(const char* name, const char* type, int w, int h, int d, float maxLen, cCellDelta3 deltas[])
    {
        // Write as stack of images
        uint8_t* dist8 = new uint8_t[w * h];

        for (int slice = 0; slice < d; slice++)
        {
            cCellDelta3* sliceDeltas = deltas + w * h * slice;

            for (int i = 0, n = w * h; i < n; i++)
                dist8[i] = EncodeU8(Dist(sliceDeltas[i]) / maxLen);

            WriteImage(name, type, w, h, dist8, slice);
        }

        delete[] dist8;
    }

#if USE_MESH
    uint32_t* BitMaskFromObjFile(int w, int h, int d, const char* sourceFile, FILE* )
    {
        FILE* meshFile = fopen(sourceFile, "r");

        if (!meshFile)
        {
            fprintf(stderr, "Can't open %s\n", sourceFile);
            return 0;
        }

        cMesh mesh;
        ReadObjFile(meshFile, &mesh);

        Bounds3f bbox = FindBounds(mesh);

        int maskSize = MaskSize(w, h, d);
        uint32_t* mask = new uint32_t[maskSize]();

        CreateBitMaskFromTriangles
        (
            (int) mesh.mPositionIndices.size() / 3,
            mesh.mPositionIndices.data(),
            mesh.mPositions.data(),
            bbox,
            w, h, d,
            mask
        );

        return mask;
    }
#endif

    void FindDeltas(int method, int w, int h, int d, cCellDelta3 deltas[], FILE* s, int cw = 1)
    {
        switch (method)
        {
        case kMethodDanielsson:
            if (s)
                fprintf(s, "\nDanielsson:\n");
            Danielsson(w, h, d, deltas, cw);
            break;
        case kMethodFastSweep:
            if (s)
                fprintf(s, "\nFast Sweep:\n");
            FastSweep(w, h, d, deltas, cw);
            break;
        case kMethodJumpFlood:
            if (s)
                fprintf(s, "\nJump Flood:\n");
            JumpFlood(w, h, d, deltas, cw);
            break;
        case kMethodBruteForce:
            if (s)
                fprintf(s, "\nBrute Force:\n");
            BruteForce(w, h, d, deltas, cw);
            break;
        default:
            fprintf(stderr, "\nUnknown method for this technique\n");
            break;
        }
    }
    
    bool GenerateDistance(const char* name, const int w, const int h, const int d, float maxLen, const uint32_t mask[], int method, FILE* s, bool check)
    {
        cCellDelta3* deltas = new cCellDelta3[w * h * d];
        
        InitDeltasFromBitMask(w, h, d, mask, deltas);

        if (s)
            fprintf(s, "Generating exterior distances\n");

        FindDeltas(method, w, h, d, deltas, s);
        
        if (s)
            PrintDeltas(w, h, d, deltas, s);

        if (check)
        {
            if (!CheckDeltas(w, h, d, deltas, stderr))
            {
                delete[] deltas;
                return false;
            }

            if (method != kMethodBruteForce)
            {
                fprintf(s, "Optimal:\n");

                cCellDelta3* optDeltas = new cCellDelta3[w * h * d];
                InitDeltasFromBitMask(w, h, d, mask, optDeltas);
                BruteForce(w, h, d, optDeltas);

                if (s)
                {
                    PrintDeltas(w, h, d, optDeltas, s);
                    fprintf(s, "Diff from optimal:\n");
                    PrintDiff(w, h, d, deltas, optDeltas, s);
                }

                PrintError(w, h, d, deltas, optDeltas, stderr);

                delete[] optDeltas;
            }
        }

        WriteDistances(name, "df", w, h, d, maxLen, deltas);
        delete[] deltas;

        return true;
    }


    ////////////////////////////////////////////////////////////////////////////
    // MARK: - Signed Distance generation
    ////////////////////////////////////////////////////////////////////////////

    enum tOutputType
    {
        kImageStandard,
        kImageDeltas,
        kImageRG,
    };

    void WriteDistances(const char* name, int w, int h, float maxLen, const cCellDelta2 deltas0[], const cCellDelta2 deltas1[], int type)
    {
        if (type == kImageStandard)
        {
            uint8_t* dist8 = new uint8_t[w * h];

            for (int i = 0, n = w * h; i < n; i++)
                dist8[i] = EncodeU8Signed((Dist(deltas0[i]) - Dist(deltas1[i])) / maxLen);

            WriteImage(name, "sdf", w, h, dist8);
            delete[] dist8;
        }
        else
        {
            uint8_t (*dist32)[4] = new uint8_t[w * h][4];

            for (int i = 0, n = w * h; i < n; i++)
            {
                if (type == kImageDeltas)
                {
                    dist32[i][0] = EncodeU8Signed((deltas0[i].x + deltas1[i].x) / maxLen);
                    dist32[i][1] = EncodeU8Signed((deltas0[i].y + deltas1[i].y) / maxLen);
                }
                else
                {
                    dist32[i][0] = EncodeU8(Dist(deltas0[i]) / maxLen);
                    dist32[i][1] = EncodeU8(Dist(deltas1[i]) / maxLen);
                }
                dist32[i][2] = 0;
                dist32[i][3] = 255;
            }

            WriteImage(name, type == kImageDeltas ? "sdf-delta" : "sdf-rg", w, h, dist32);
            delete[] dist32;
        }
    }

    bool GenerateSignedDistance(const char* name, int w, int h, float maxLen, uint32_t mask[], int method, int imageType, FILE* s, bool check)
    {
        cCellDelta2* deltas0 = new cCellDelta2[w * h];
        cCellDelta2* deltas1 = new cCellDelta2[w * h];

        InitDeltasFromBitMask(w, h, mask, deltas0, false);
        InitDeltasFromBitMask(w, h, mask, deltas1, true);

        if (s)
        {
            PrintDeltas(w, h, deltas0, s);
            PrintDeltas(w, h, deltas1, s);

            fprintf(s, "Generating exterior and interior distances\n");
        }

        FindDeltas(method, w, h, deltas0, s);
        FindDeltas(method, w, h, deltas1, s);

        if (s)
        {
            PrintDeltas(w, h, deltas0, s);
            PrintDeltas(w, h, deltas1, s);
        }

        if (check)
        {
            if (!CheckDeltas(w, h, deltas0, stderr))
                return false;

            if (!CheckDeltas(w, h, deltas1, stderr))
                return false;

            // Not much point checking vs. optimal as covered by distance generation path
        }

        WriteDistances(name, w, h, maxLen, deltas0, deltas1, imageType);

        delete[] deltas0;
        delete[] deltas1;

        return true;
    }

    void WriteDistances(const char* name, int w, int h, int d, float maxLen, cCellDelta3 deltas0[], cCellDelta3 deltas1[], int type)
    {
        if (type == kImageStandard)
        {
            uint8_t* dist8 = new uint8_t[w * h];

            for (int slice = 0; slice < d; slice++)
            {
                cCellDelta3* sliceDeltas0 = deltas0 + w * h * slice;
                cCellDelta3* sliceDeltas1 = deltas1 + w * h * slice;

                for (int i = 0, n = w * h; i < n; i++)
                    dist8[i] = EncodeU8Signed((Dist(sliceDeltas0[i]) - Dist(sliceDeltas1[i])) / maxLen);

                WriteImage(name, "sdf", w, h, dist8, slice);
            }

            delete[] dist8;
        }
        else
        {
            uint8_t (*dist32)[4] = new uint8_t[w * h][4];

            for (int slice = 0; slice < d; slice++)
            {
                cCellDelta3* sliceDeltas0 = deltas0 + w * h * slice;
                cCellDelta3* sliceDeltas1 = deltas1 + w * h * slice;

                for (int i = 0, n = w * h; i < n; i++)
                {
                    if (type == kImageDeltas)
                    {
                        dist32[i][0] = EncodeU8Signed((sliceDeltas0[i].x + sliceDeltas1[i].x) / maxLen);
                        dist32[i][1] = EncodeU8Signed((sliceDeltas0[i].y + sliceDeltas1[i].y) / maxLen);
                        dist32[i][2] = EncodeU8Signed((sliceDeltas0[i].z + sliceDeltas1[i].z) / maxLen);
                    }
                    else
                    {
                        dist32[i][0] = EncodeU8(Dist(sliceDeltas0[i]) / maxLen);
                        dist32[i][1] = EncodeU8(Dist(sliceDeltas1[i]) / maxLen);
                        dist32[i][2] = 0;
                    }
                    dist32[i][3] = 255;
                }

                WriteImage(name, type == kImageDeltas ? "sdf-delta" : "sdf-rg", w, h, dist32, slice);
            }

            delete[] dist32;
        }
    }

    bool GenerateSignedDistance(const char* name, const int w, const int h, const int d, float maxLen, const uint32_t mask[], int method, int imageType, FILE* s, bool check)
    {
        cCellDelta3* deltas0 = new cCellDelta3[w * h * d];
        cCellDelta3* deltas1 = new cCellDelta3[w * h * d];

        InitDeltasFromBitMask(w, h, d, mask, deltas0, false);
        InitDeltasFromBitMask(w, h, d, mask, deltas1, true);

        if (s)
            fprintf(s, "Generating exterior and interior distances\n");

        FindDeltas(method, w, h, d, deltas0, s);
        FindDeltas(method, w, h, d, deltas1, s);

        if (s)
        {
            PrintDeltas(w, h, d, deltas0, s);
            PrintDeltas(w, h, d, deltas1, s);
        }

        if (check)
        {
            if (!CheckDeltas(w, h, d, deltas0, stderr))
                return false;

            if (!CheckDeltas(w, h, d, deltas1, stderr))
                return false;
        }

        WriteDistances(name, w, h, d, maxLen, deltas0, deltas1, imageType);

        delete[] deltas0;
        delete[] deltas1;

        return true;
    }



    ////////////////////////////////////////////////////////////////////////////
    // MARK: - Signed Distance generation using border technique
    ////////////////////////////////////////////////////////////////////////////

    void WriteDistances(const char* name, int w, int h, float maxLen, const uint32_t mask[], const cCellDelta2 deltas[], int type)
    {
        if (type == kImageStandard)
        {
            uint8_t* dist8 = new uint8_t[w * h];

            for (int y = 0; y < h; y++)
            {
                const cCellDelta2* deltasRow = deltas + w * y;
                uint8_t*           dist8Row  = dist8  + w * y;

                for (int j = 0; j < w; j += 32)
                {
                    uint32_t m = *mask++;

                    int n = 32;
                    if (w - j < n)
                        n = w - j;

                    for (int i = 0; i < n; i++)
                    {
                        float d = Dist(deltasRow[j + i]) / maxLen;
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
                const cCellDelta2* deltasRow     = deltas + w * y;
                uint8_t          (*dist32Row)[4] = dist32 + w * y;

                if (type == kImageDeltas)
                {
                    for (int x = 0; x < w; x++)
                    {
                        dist32Row[x][0] = EncodeU8Signed(deltasRow[x].x / maxLen);
                        dist32Row[x][1] = EncodeU8Signed(deltasRow[x].y / maxLen);
                        dist32Row[x][2] = 0;
                        dist32Row[x][3] = 255;
                    }
                }
                else
                    for (int j = 0; j < w; j += 32)
                    {
                        uint32_t m = *mask++;

                        int n = 32;
                        if (w - j < n)
                            n = w - j;

                        for (int i = 0; i < n; i++)
                        {
                            dist32Row[j + i][ m & 1     ] = EncodeU8(Dist(deltasRow[j + i]) / maxLen);
                            dist32Row[j + i][(m & 1) ^ 1] = 0;
                            dist32Row[j + i][2] = 0;
                            dist32Row[j + i][3] = 255;

                            m >>= 1;
                        }
                    }
            }

            WriteImage(name, type == kImageDeltas ? "sdf-delta" : "sdf-rg", w, h, dist32);
            delete[] dist32;
        }
    }

    bool GenerateSignedDistanceBorder(const char* name, int w, int h, float maxLen, uint32_t mask[], int method, int imageType, FILE* s, bool check)
    {
        cCellDelta2* deltas = new cCellDelta2[w * h];

        InitBorderDeltasFromBitMask(w, h, mask, deltas);

        if (s)
            PrintDeltas(w, h, deltas, s);

        if (check && !CheckBorderDeltas(w, h, mask, deltas, stderr))
            return false;

        if (s)
            fprintf(s, "Generating signed distances from border init\n");

        FindDeltas(method, w, h, deltas, s, 2);
        
        if (s)
            PrintDeltas(w, h, deltas, s);

        if (check)
        {
            if (!CheckDeltas(w, h, deltas, stderr, 2))
                return false;

            // Not much point checking vs. optimal as covered by distance generation path
        }

        WriteDistances(name, w, h, maxLen * 2, mask, deltas, imageType);
        delete[] deltas;

        return true;
    }

    void WriteDistances(const char* name, int w, int h, int d, float maxLen, const uint32_t mask[], cCellDelta3 deltas[], int type)
    {
        if (type == kImageStandard)
        {
            // Write as images
            uint8_t* dist8 = new uint8_t[w * h];

            for (int slice = 0; slice < d; slice++)
            {
                for (int y = 0; y < h; y++)
                {
                    const cCellDelta3* deltasRow = deltas + w * y;
                    uint8_t*            dist8Row = dist8  + w * y;

                    for (int j = 0; j < w; j += 32)
                    {
                        uint32_t m = *mask++;

                        int n = 32;
                        if (w - j < n)
                            n = w - j;

                        for (int i = 0; i < n; i++)
                        {
                            float d = Dist(deltasRow[j + i]) / maxLen;
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

            for (int slice = 0; slice < d; slice++)
            {
                for (int y = 0; y < h; y++)
                {
                    const cCellDelta3* deltasRow     = deltas + w * y;
                    uint8_t          (*dist32Row)[4] = dist32 + w * y;

                    if (type == kImageDeltas)
                    {
                        for (int x = 0; x < w; x++)
                        {
                            dist32Row[x][0] = EncodeU8Signed(deltasRow[x].x / maxLen);
                            dist32Row[x][1] = EncodeU8Signed(deltasRow[x].y / maxLen);
                            dist32Row[x][2] = EncodeU8Signed(deltasRow[x].z / maxLen);
                            dist32Row[x][3] = 255;
                        }
                    }
                    else
                    {
                        for (int j = 0; j < w; j += 32)
                        {
                            uint32_t m = *mask++;

                            int n = 32;
                            if (w - j < n)
                                n = w - j;

                            for (int i = 0; i < n; i++)
                            {
                                dist32Row[j + i][ m & 1     ] = EncodeU8(Dist(deltasRow[j + i]) / maxLen);
                                dist32Row[j + i][(m & 1) ^ 1] = 0;
                                dist32Row[j + i][2] = 0;
                                dist32Row[j + i][3] = 255;

                                m >>= 1;
                            }
                        }
                    }
                }

                WriteImage(name, type == kImageDeltas ? "sdf-delta" : "sdf-rg", w, h, dist32, slice);
                deltas += w * h;
            }

            delete[] dist32;
        }
    }

    bool GenerateSignedDistanceBorder(const char* name, const int w, const int h, const int d, float maxLen, const uint32_t mask[], int method, int imageType, FILE* s, bool check)
    {
        cCellDelta3* deltas = new cCellDelta3[w * h * d];

        InitBorderDeltasFromBitMask(w, h, d, mask, deltas);

        if (s)
            PrintDeltas(w, h, d, deltas, s);

        if (check && !CheckBorderDeltas(w, h, d, mask, deltas, stderr))
            return false;

        if (s)
            fprintf(s, "Generating signed distances from border init\n");

        FindDeltas(method, w, h, d, deltas, s, 2);

        if (s)
            PrintDeltas(w, h, d, deltas, s);

        if (check && !CheckDeltas(w, h, d, deltas, stderr, 2))
            return false;

        WriteDistances(name, w, h, d, maxLen * 2, mask, deltas, imageType);
        delete[] deltas;

        return true;
    }




    ////////////////////////////////////////////////////////////////////////////
    // MARK: - Occlusion testing
    ////////////////////////////////////////////////////////////////////////////

    void GenerateOcclusion(const char* name, const int w, const int h, float dirW[], int method, FILE* s)
    {
        OcclusionSweep(w, h, dirW);

        const int imageSize = w * h;
        float* dirW4[4] = { dirW, dirW + imageSize, dirW + imageSize * 2, dirW + imageSize * 3 };

        if (s)
            for (int i = 0; i < 4; i++)
            {
                fprintf(s, "quadrant %d:\n", i);
                PrintDistances(w, h, dirW4[i], s);
            }

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

            if (method == 2)
                aoDir[i][2] = aoDir[i][1];
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

    void GenerateOcclusion(const char* name, const int w, const int h, const int d, float* dirW, int method, FILE* s)
    {
        const size_t volumeStride = w * h * d;

        if (s)
            fprintf(s, "\nMethod %d:\n", method);

        if (method & 1)
        {
            float* outDirW = new float[8 * volumeStride];

            OcclusionSweep(w, h, d, dirW, outDirW);
            memcpy(dirW, outDirW, 8 * volumeStride * sizeof(float));

            delete[] outDirW;
        }
        else
            OcclusionSweep(w, h, d, dirW);

        if (s)
            for (int i = 0; i < 8; i++)
            {
                fprintf(s, "octant %d:\n", i);
                PrintDistances(w, h, d, dirW + volumeStride * i, s);
            }

        // Write as images
        const size_t sliceStride = w * h;
        uint8_t* ao         = new uint8_t[sliceStride];
        uint8_t (*aoDir)[4] = new uint8_t[sliceStride][4];
        
        for (int slice = 0; slice < d; slice++)
        {
            float* sDirW[8];

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

                if (method == 2)
                    aoDir[i][3] = 255;
            }
            
            WriteImage(name, "dirW", w, h, aoDir, slice);
        }

        delete[] ao;
        delete[] aoDir;
    }
    
#if USE_MESH
    float* DirWFromFile(const int w, const int h, const int d, const char* sourceFile)
    {
        FILE* meshFile = fopen(sourceFile, "r");

        if (!meshFile)
        {
            fprintf(stderr, "Can't open %s\n", sourceFile);
            return 0;
        }

        cMesh mesh;
        ReadObjFile(meshFile, &mesh);

        const size_t volumeStride = w * h * d;

        float* dirW = new float[8 * volumeStride]();

        Bounds3f bbox = FindBounds(mesh);
        bbox = FindAOBounds(0.5f, bbox);

        CreateDirW8FromTriangles
        (
            (int) mesh.mPositionIndices.size() / 3,
            mesh.mPositionIndices.data(),
            mesh.mPositions.data(),
            bbox,
            w, h, d,
            dirW
        );

        return dirW;
    }
#endif
}

////////////////////////////////////////////////////////////////////////////
// MARK: - Tool
////////////////////////////////////////////////////////////////////////////

namespace
{
    inline const char* Next(int& argc, const char**& argv)
    {
        const char* result = argc > 0 ? argv[0] : 0;
        argc--; argv++;
        return result;
    }
    inline bool NextIsOption(int argc, const char* argv[])
    {
        return argc > 0 && argv[0][0] == '-' && isalpha(argv[0][1]);
    }
    inline bool NextIsArg(int argc, const char* argv[])
    {
        return argc > 0 && (argv[0][0] != '-' || !isalpha(argv[0][1]));
    }

    void GetFileName(char* buffer, size_t bufferSize, const char* path)
    {
        const char* lastSlash = strrchr(path, '/');
        if (!lastSlash)
            lastSlash = strrchr(path, '\\');

        if (lastSlash)
            strlcpy(buffer, lastSlash + 1, bufferSize);
        else
            strlcpy(buffer, path, bufferSize);

        char* lastDot = strrchr(buffer, '.');

        if (lastDot)
            *lastDot = 0;
    }
}

int main(int argc, const char* argv[])
{
    if (argc == 1)
    {
        printf
        (
            "Distance/Occlusion generation tool\n"
            "\n"
            "Options:\n"
            "  -d: generate distance field. Methods: 0, Danielsson; 1, Fast Sweep; 2, Jump Flood; 3, Brute Force; 4, Chamfer; 5, Felzenszwalb\n"
            "  -s: generate signed distance field in a single pass. Methods: as for -d but Chamfer/Felzenszwalb unsupported.\n"
            "  -S: generate signed distance field via orthodox two distance field passes.\n"
            "  -x: generate occlusion. Methods 0: standard, 1: no self-occlusion, 2: directional components only\n"
            "\n"
            "  -m <int>: select method variant for above\n"
            "\n"
            "  -w <w:int> [<h:int>]        : set dimensions of image source\n"
            "  -W <w:int> [<h:int> <d:int>]: set dimensions of volume source\n"
            "  -p [<count:int> <seed:int>] : add random points to image/volume (default)\n"
            "  -b <sides:int>              : add test box with the given number of sides to image/volume\n"
        #if USE_STB
            "\n"
            "  -f <path>: specify input image\n"
            "  -t <int> : specify 0-255 threshold for creating mask from image\n"
            "  -r       : reverse so white is occupied rather than black\n"
        #endif
        #if USE_MESH
            "\n"
            "  -F <path>: use given obj file to define mask.\n"
        #endif
            "\n"
            "  -o <name>: set base name for output file(s)\n"
        #if USE_STB || USE_MESH
            "  -v       : log detailed output\n"
        #endif
            "  -c       : run checks on output\n"
        );
        
        return -1;
    }

    Next(argc, argv); // command name

    enum tAlgorithm
    {
        kNone,
        kDistance,
        kSignedDistance,
        kSignedDistanceBorder,
        kOcclusion,
        kMaxAlgorithms
    };

    int algorithm = kNone;
    int method = 0;
    int imageType = kImageStandard;
    int dim[3] = { 0, 0, 0 };
    float maxLen = 0.0f;
    bool volume = false;

    int pointsSeed = 12345;
    int numPoints = 0;
    int numSides = 0;

    char outName[1024] = "out";
    bool haveName = false;
    int check = 0;

#if USE_STB || USE_MESH
    int threshold = 128;
    bool lessThan = true;
    const char* sourceFile = 0;
    FILE* out = 0;
#else
    FILE* out = stdout;
#endif

    while (NextIsOption(argc, argv))
    {
        const char* option = Next(argc, argv) + 1;

        switch (option[0])
        {
        case 'd':
            algorithm = kDistance;
            break;
        case 's':
            algorithm = kSignedDistanceBorder;
            break;
        case 'S':
            algorithm = kSignedDistance;
            break;
        case 'x':
            algorithm = kOcclusion;
            break;
        case 'v':
            out = stdout;
            break;
        case 'c':
            check = 1;
            break;
        case 'm':
            if (NextIsArg(argc, argv))
                method = atoi(Next(argc, argv));
            break;
        case 'i':
            if (NextIsArg(argc, argv))
            {
                imageType = atoi(Next(argc, argv));
                
                if (maxLen == 0 && imageType == kImageDeltas)
                    maxLen = 127;   // default to 1:1
            }
            break;
        case 'w':
        case 'W':
            for (int i = 0; i < 3 && NextIsArg(argc, argv); i++)
                dim[i] = atoi(Next(argc, argv));

            volume = (option[0] == 'W');
            break;
        case 'p':
            if (NextIsArg(argc, argv))
                numPoints = atoi(Next(argc, argv));
            else
                numPoints = (dim[0] + dim[1] + dim[2]) / 8;

            if (NextIsArg(argc, argv))
                pointsSeed = atoi(Next(argc, argv));
            break;
        case 'b':
            if (NextIsArg(argc, argv))
                numSides = atoi(Next(argc, argv));
            else
                numSides = 2;
            break;
    #if USE_STB || USE_MESH
        case 't':
            if (NextIsArg(argc, argv))
                threshold = atoi(Next(argc, argv));
            break;
        case 'r':
            lessThan = !lessThan;
            break;
        case 'l':
            if (NextIsArg(argc, argv))
                maxLen = atoi(Next(argc, argv));
            break;
        case 'f':
        case 'F':
            if (NextIsArg(argc, argv))
                sourceFile = Next(argc, argv);

            volume = volume || (option[0] == 'F');
            break;
    #endif
        case 'o':
            if (NextIsArg(argc, argv))
            {
                strncpy(outName, Next(argc, argv), sizeof(outName));
                haveName = true;
            }
            break;
        default:
            fprintf(stderr, "Unknown option: -%s\n", option);
            return -1;
        }
    }

    if (argc > 0)
    {
        fprintf(stderr, "Unknown argument: %s\n", argv[0]);
        return -1;
    }

    if (dim[0] == 0)
        dim[0] = volume ? 32 : 256;

    if (dim[1] == 0)
        dim[1] = dim[0];

    if (volume && dim[2] == 0)
        dim[2] = dim[1];

    if (!sourceFile && numSides == 0 && numPoints == 0)
        numPoints = (dim[0] + dim[1] + dim[2]) / 8;

    uint32_t* mask = 0;
    float* dirW = 0;

    if (out)
        fprintf(out, "Generating presence mask\n");

    if (!volume)
    {
    #if USE_STB
        if (sourceFile)
            mask = BitMaskFromImageFile(sourceFile, threshold, lessThan, dim[0], dim[1]);
        else
    #endif
            mask = CreateBitMask(dim[0], dim[1]);

        if (!mask)
            return -1;

        if (numSides > 0)
            BitMaskAddBlock(dim[0], dim[1], numSides, mask);
        if (numPoints > 0)
            BitMaskAddPoints(dim[0], dim[1], numPoints, pointsSeed, mask);

        if (out)
            PrintMask(dim[0], dim[1], mask, out);

        if (algorithm == kOcclusion)
        {
            dirW = new float[4 * dim[0] * dim[1]];
            InitDirWFromBitMask(dim[0], dim[1], mask, dirW);
            delete[] mask;
            mask = 0;
        }
    }
    else
    {
    #if USE_MESH
        if (sourceFile)
        {
            if (algorithm == kOcclusion)
                dirW = DirWFromFile(dim[0], dim[1], dim[2], sourceFile);    // construct DirW directly so we can take advantage of triangle normal info
            else
                mask = BitMaskFromObjFile(dim[0], dim[1], dim[2], sourceFile, out);
        }
        else
    #endif
            mask = CreateBitMask(dim[0], dim[1], dim[2]);

        if (!mask && !dirW)
            return -1;
        
        if (mask)
        {
            if (numSides > 0)
                BitMaskAddBlock(dim[0], dim[1], dim[2], numSides, mask);
            if (numPoints > 0)
                BitMaskAddPoints(dim[0], dim[1], dim[2], numPoints, pointsSeed, mask);
            if (out)
                PrintMask(dim[0], dim[1], dim[1], mask, out);

            if (algorithm == kOcclusion)
            {
                dirW = new float[8 * dim[0] * dim[1] * dim[2]];
                InitDirWFromBitMask(dim[0], dim[1], dim[2], mask, dirW);
                delete[] mask;
                mask = 0;
            }
        }
    }

    if (maxLen == 0.0f)
        maxLen = 0.25f * (dim[0] >= dim[1] ? dim[0] : dim[1]);

    bool success = true;

    if (out)
        fprintf(out, "Generating results\n");

#if USE_STB || USE_MESH
    if (!haveName && sourceFile)
        GetFileName(outName, sizeof(outName), sourceFile);
#endif

    switch (algorithm)
    {
    case kDistance:
        if (dim[2] == 0)
            success = GenerateDistance(outName, dim[0], dim[1], maxLen, mask, method, out, check);
        else
            success = GenerateDistance(outName, dim[0], dim[1], dim[2], maxLen, mask, method, out, check);
        break;

    case kSignedDistance:
        if (dim[2] == 0)
            success = GenerateSignedDistance(outName, dim[0], dim[1], maxLen, mask, method, imageType, out, check);
        else
            success = GenerateSignedDistance(outName, dim[0], dim[1], dim[2], maxLen, mask, method, imageType, out, check);
        break;

    case kSignedDistanceBorder:
        if (dim[2] == 0)
            success = GenerateSignedDistanceBorder(outName, dim[0], dim[1], maxLen, mask, method, imageType, out, check);
        else
            success = GenerateSignedDistanceBorder(outName, dim[0], dim[1], dim[2], maxLen, mask, method, imageType, out, check);
        break;

    case kOcclusion:
        if (dim[2] == 0)
            GenerateOcclusion(outName, dim[0], dim[1], dirW, method, out);
        else
            GenerateOcclusion(outName, dim[0], dim[1], dim[2], dirW, method, out);
        break;

    default:
        fprintf(stderr, "Please specify an algorithm (-d/-s/-S/-x)\n");
        success = false;
    }

    delete[] mask;
    delete[] dirW;

    return success ? 0 : -1;
}

