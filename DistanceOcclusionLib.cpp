//------------------------------------------------------------------------------
// Purpose: Distance Field and Occlusion Generation
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#include "DistanceOcclusionLib.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace DOL;

// Declarations


namespace
{
    // To get the full INT16_MAX range, Closer() needs to
    // use int64_t for its partial results in the Closer(x, y, z) case.
    // We just limit the max to avoid this.
    const int16_t kMaxDelta2 = INT16_MAX;
    const int16_t kMaxDelta3 = 25000;

    inline uint32_t NextSeed(uint32_t seed)
    {
        return seed * 0x278dde6d;
    }
}


// --- 2D Images ---------------------------------------------------------------

void DOL::CreateBitMask(int w, int h, uint32_t seed, uint32_t mask[])
{
    int scale = (seed >> 4);
    seed &= 15;
    
    int count = (w + h) * (1 + scale) / 8;
    int ws = MaskSize(w);
    
    for (int i = 0; i < count; i++)
    {
        seed = NextSeed(seed);
        int r = seed & ~0x80000000;
        
        int s = r / w;
        int t = s / h;
        
        int x = r - w * s;
        int y = s - h * t;
        
        assert(x < w);
        assert(y < h);
        
        int fx = x & 31;
        x >>= 5;
        
        int index = y * ws + x;
        mask[index] |= (1 << fx);
    }
}


void DOL::InitDFFromBitMask(int w, int h, const uint32_t mask[], float distances[])
{
    // mask is expected to have rows padded to whole numbers of uint32_ts
    int s = w * h;
    for (int i = 0; i < s; i++)
        distances[i] = FLT_MAX;
    
    int maskStride = MaskSize(w);
    
    for (int y = 0; y < h; y++)
        for (int j = 0; j < maskStride; j++)
        {
            uint32_t m = *mask++;
            
            if (m == 0)
                continue;
            
            float* block = distances + y * w + 32 * j;
            
            for (int i = 0; i < 32; i++)
            {
                if (m & 1)
                    block[i] = 0.0f;
                
                m >>= 1;
            }
        }
}

void DOL::InitDFFromBitMask(int w, int h, const uint32_t mask[], int32_t distances[])
{
    // mask is expected to have rows padded to whole numbers of uint32_ts
    int s = w * h;
    for (int i = 0; i < s; i++)
        distances[i] = INT32_MAX;

    int maskStride = MaskSize(w);
    
    for (int y = 0; y < h; y++)
    {
        for (int j = 0; j < maskStride; j++)
        {
            uint32_t m = *mask++;
                        
            if (m == 0)
                continue;
            
            int32_t* block = distances + y * w + 32 * j;
            
            for (int i = 0; i < 32; i++)
            {
                if (m & 1)
                    block[i] = 0;
                
                m >>= 1;
            }
        }
    }
}

void DOL::Chamfer(int w, int h, float distances[])
{
    // according to Borgefors, using these values leads to more
    // stable and accurate results than 1, sqrt(2) with double precision even
    // ceil(d) is still a metric, 
    // G. Borgefors.
    // Distance transformations in digital images.
    
    const float d0 = 3;
    const float d1 = 4;
    
    float* row0 = distances;
    
    for (int y = 1; y < h; y++)
    {
        float* row1 = row0 + w;
        
        if (row1[0] > row0[0] + d0)
            row1[0] = row0[0] + d0;
        if (row1[0] > row0[1] + d1)
            row1[0] = row0[1] + d1;
        
        for (int x = 1; x < w - 1; x++)
        {
            if (row1[x] > row0[x - 1] + d1)
                row1[x] = row0[x - 1] + d1;
            if (row1[x] > row0[x    ] + d0)
                row1[x] = row0[x    ] + d0;
            if (row1[x] > row0[x + 1] + d1)
                row1[x] = row0[x + 1] + d1;
            
            if (row1[x] > row1[x - 1] + d0)
                row1[x] = row1[x - 1] + d0;
        }
        
        if (row1[w - 1] > row0[w - 2] + d1)
            row1[w - 1] = row0[w - 2] + d1;
        if (row1[w - 1] > row0[w - 1] + d0)
            row1[w - 1] = row0[w - 1] + d0;
        
        if (row1[w - 1] > row1[w - 2] + d0)
            row1[w - 1] = row1[w - 2] + d0;
        
        row0 = row1;
    }

    row0 = distances + h * w - w;
    
    for (int y = h - 2; y >= 0; y--)
    {
        float* row1 = row0 - w;
        
        if (row1[w - 1] > row0[w - 1] + d0)
            row1[w - 1] = row0[w - 1] + d0;
        if (row1[w - 1] > row0[w - 2] + d1)
            row1[w - 1] = row0[w - 2] + d1;
        
        for (int x = w - 2; x > 0; x--)
        {
            if (row1[x] > row0[x + 1] + d1)
                row1[x] = row0[x + 1] + d1;
            if (row1[x] > row0[x    ] + d0)
                row1[x] = row0[x    ] + d0;
            if (row1[x] > row0[x - 1] + d1)
                row1[x] = row0[x - 1] + d1;
            
            if (row1[x] > row1[x + 1] + d0)
                row1[x] = row1[x + 1] + d0;
        }
        
        if (row1[0] > row0[1] + d1)
            row1[0] = row0[1] + d1;
        if (row1[0] > row0[0] + d0)
            row1[0] = row0[0] + d0;
        
        if (row1[0] > row1[1] + d0)
            row1[0] = row1[1] + d0;
        
        row0 = row1;
    }
}

void DOL::Chamfer(int w, int h, int32_t distances[])
{
    // according to Borgefors, using these values leads to more
    // stable and accurate results than 1, sqrt(2) with double precision even
    // ceil(d) is still a metric, 
    // G. Borgefors.
    // Distance transformations in digital images.
    
    const int d0 = 3;
    const int d1 = 4;
    
    int32_t* row0 = distances;
    
    for (int y = 1; y < h; y++)
    {
        int32_t* row1 = row0 + w;
        
        if (row1[0] > row0[0] + d0)
            row1[0] = row0[0] + d0;
        if (row1[0] > row0[1] + d1)
            row1[0] = row0[1] + d1;
        
        for (int x = 1; x < w - 1; x++)
        {
            if (row1[x] > row0[x - 1] + d1)
                row1[x] = row0[x - 1] + d1;
            if (row1[x] > row0[x    ] + d0)
                row1[x] = row0[x    ] + d0;
            if (row1[x] > row0[x + 1] + d1)
                row1[x] = row0[x + 1] + d1;
            
            if (row1[x] > row1[x - 1] + d0)
                row1[x] = row1[x - 1] + d0;
        }
        
        if (row1[w - 1] > row0[w - 2] + d1)
            row1[w - 1] = row0[w - 2] + d1;
        if (row1[w - 1] > row0[w - 1] + d0)
            row1[w - 1] = row0[w - 1] + d0;
        
        if (row1[w - 1] > row1[w - 2] + d0)
            row1[w - 1] = row1[w - 2] + d0;
        
        row0 = row1;
    }
    
    row0 = distances + h * w - w;
    
    for (int y = h - 2; y >= 0; y--)
    {
        int32_t* row1 = row0 - w;
        
        if (row1[w - 1] > row0[w - 1] + d0)
            row1[w - 1] = row0[w - 1] + d0;
        if (row1[w - 1] > row0[w - 2] + d1)
            row1[w - 1] = row0[w - 2] + d1;
        
        for (int x = w - 2; x > 0; x--)
        {
            if (row1[x] > row0[x + 1] + d1)
                row1[x] = row0[x + 1] + d1;
            if (row1[x] > row0[x    ] + d0)
                row1[x] = row0[x    ] + d0;
            if (row1[x] > row0[x - 1] + d1)
                row1[x] = row0[x - 1] + d1;
            
            if (row1[x] > row1[x + 1] + d0)
                row1[x] = row1[x + 1] + d0;
        }
        
        if (row1[0] > row0[1] + d1)
            row1[0] = row0[1] + d1;
        if (row1[0] > row0[0] + d0)
            row1[0] = row0[0] + d0;
        
        if (row1[0] > row1[1] + d0)
            row1[0] = row1[1] + d0;
        
        row0 = row1;
    }
}

// To adapt to 3D: need new chamfer approx values. See DistanceTransform.pdf
// To adapt to height field: ...

void DOL::InitDeltaFromBitMask(int w, int h, const uint32_t mask[], cCellDelta2 delta[])
{
    cCellDelta2 maxDelta = { kMaxDelta2, kMaxDelta2 };
    cCellDelta2 minDelta = { 0, 0 };
    
    int s = w * h;
    for (int i = 0; i < s; i++)
        delta[i] = maxDelta;
    
    const int maskStride = MaskSize(w);
    
    for (int y = 0; y < h; y++)
        for (int j = 0; j < maskStride * 32; j += 32)
        {
            uint32_t m = (*mask++);
            
            if (m == 0)
                continue;
            
            cCellDelta2* block = delta + y * w + j;
            
            for (int i = 0; i < 32; i++)
            {
                if (m & 1)
                    block[i] = minDelta;
                
                m >>= 1;
            }
        }
}


// Theory: Danielsson sweep is the same as the fast sweeping except
// it wraps the alternating x scans into one vertical sweep.
// Despite claims in papers, it seems its mem access density is approx 6wh,
// vs fs' 4 x 2wh = 8wh.
// in 3D, Danielsson is 8wh * d + wh x d-1 x 2 =~ 10whd
// vs 8 x 3 x whd(!)

#define ONE_WAY_X 0
// this doesn't make a difference, because first rows are dealt with by
// the last row of the alternate sweep
#define DO_FIRST_ROWS 1

namespace
{
    /// Returns true if d0 is closer than d1 
    inline bool Closer(cCellDelta2 d0, cCellDelta2 d1)
    {
        int dw0 = d0.x * d0.x + d0.y * d0.y;
        int dw1 = d1.x * d1.x + d1.y * d1.y;
        
        return dw0 > dw1;
    }
}

void DOL::Danielsson(int w, int h, cCellDelta2 delta[])
{
    // store dx^2, dy^2 to nearest point. distance = sqrt(dx^2 + dy^2).
    // - mostly accurate, doesn't overestimate like Chamfer.
    // Cui99 CM99 fix errors
    // O. Cuisenaire and B. Macq.
    // Fast euclidean distance transformations by propagation using multiple neighbourhoods.    

    cCellDelta2* row0 = delta;
    cCellDelta2 cell;

#if DO_FIRST_ROWS
    for (int i = 1; i < w; i++)
    {
        cell.x = row0[i - 1].x - 1;
        cell.y = row0[i - 1].y;

        if (Closer(row0[i], cell))
            row0[i] = cell;
    }
#endif
    
    for (int j = 1; j < h; j++)
    {
        cCellDelta2* row1 = row0 + w;
        
        for (int i = 0; i < w; i++)
        {
            cell.x = row0[i].x;
            cell.y = row0[i].y - 1;

            if (Closer(row1[i], cell))
                row1[i] = cell;
        }
        
        for (int i = 1; i < w; i++)
        {
            cell.x = row1[i - 1].x - 1;
            cell.y = row1[i - 1].y;

            if (Closer(row1[i], cell))
                row1[i] = cell;
        }
#if !ONE_WAY_X
        for (int i = w - 2; i >= 0; i--)
        {
            cell.x = row1[i + 1].x + 1;
            cell.y = row1[i + 1].y;
            
            if (Closer(row1[i], cell))
                row1[i] = cell;
        }
#endif
        
        row0 = row1;
    }

    row0 = delta + h * w - w;
    
#if DO_FIRST_ROWS
    for (int i = w - 2; i >= 0; i--)
    {
        cell.x = row0[i + 1].x + 1;
        cell.y = row0[i + 1].y;

        
        if (Closer(row0[i], cell))
            row0[i] = cell;
    }
#endif
    
    for (int j = h - 1; j > 0; j--)
    {
        cCellDelta2* row1 = row0 - w;
        
        for (int i = 0; i < w; i++)
        {
            cell.x = row0[i].x;
            cell.y = row0[i].y + 1;
            
            if (Closer(row1[i], cell))
                row1[i] = cell;
        }
        
#if !ONE_WAY_X
        for (int i = 1; i < w; i++)
        {
            cell.x = row1[i - 1].x - 1;
            cell.y = row1[i - 1].y;
            
            if (Closer(row1[i], cell))
                row1[i] = cell;
        }
#endif
        
        for (int i = w - 2; i >= 0; i--)
        {
            cell.x = row1[i + 1].x + 1;
            cell.y = row1[i + 1].y;
            
            if (Closer(row1[i], cell))
                row1[i] = cell;
        }
        
        row0 = row1;
    }
}

void DOL::FastSweep(int w, int h, cCellDelta2 delta[])
{
    // Fast sweep approach -- sweep each diagonal in turn.
    int lastRow = w * (h - 1);
    
    int sdx[4] = { +1, -1, +1, -1 };
    int sdy[4] = { +1, +1, -1, -1 };
    
    int ib[2] = { 0, w - 1 };
    int ie[2] = { w,    -1 };
    
    int rowStart[2] = { 0, lastRow };
    cCellDelta2 cell;

    for (int sweep = 0; sweep < 4; sweep++)
    {
        const int dx = sdx[sweep];
        const int dy = sdy[sweep];

        const int sweepX = (sweep >> 0) & 1;
        const int sweepY = (sweep >> 1) & 1;
        
        cCellDelta2* row0 = delta + rowStart[sweepY];
        
        for (int j = 1; j < h; j++)
        {
            cCellDelta2* row1 = row0 + w * dy;
        
            for (int i = ib[sweepX]; i != ie[sweepX]; i += dx)
            {
                cell.x = row0[i].x;
                cell.y = row0[i].y - dy;
                
                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }

            for (int i = ib[sweepX] + dx; i != ie[sweepX]; i += dx)
            {
                cell.x = row1[i - dx].x - dx;
                cell.y = row1[i - dx].y;

                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }
            
            row0 = row1;
        }
    }
}

// Brute force -- for checking only!
void DOL::FindOptimalDeltas(int w, int h, cCellDelta2 delta[])
{
    for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++)
        {
            cCellDelta2* cell = delta + i * w + j;
            if (cell->x != 0 || cell->y != 0)
                continue;

            for (int y = 0; y < h; y++)
                for (int x = 0; x < w; x++)
                {
                    cCellDelta2* cell = delta + y * w + x;
                    
                    int d0 = (i - y) * (i - y) + (j - x) * (j - x); // d2 to i, j
                    int d1 = cell->x * cell->x + cell->y * cell->y; // current
                    
                    if (d0 < d1)
                    {
                        cell->x = j - x;
                        cell->y = i - y;
                    }
                }
            
        }
}



namespace
{
    void FillBlock(int rs, int x0, int x1, int y0, int y1, float* dirW[4])
    {
        for (int y = y0; y <= y1; y++)
            for (int x = x0; x <= x1; x++)
            {
                int i = y * rs + x;
                
                for (int k = 0; k < 4; k++)
                    dirW[k][i] = 1.0f;
            }
    }

    uint32_t LeftMask(int n)
    {
        uint32_t result = 0;

        for ( ; n < 32; n++)
            result |= 1 << n;

        return result;
    }

    uint32_t RightMask(int n)
    {
        uint32_t result = 0;

        for (int i = 0; i <= n; i++)
            result |= 1 << i;

        return result;
    }

    void FillMaskBlock(int rs, int x0, int x1, int y0, int y1, uint32_t mask[])
    {
        int sx0 = x0 / 32;
        int sx1 = x1 / 32;

        uint32_t mask0 = LeftMask (x0 - sx0 * 32);
        uint32_t mask1 = RightMask(x1 - sx1 * 32);

        if (sx0 == sx1)
        {
            mask0 = mask0 & mask1;
            mask1 = 0;
        }

        for (int y = y0; y <= y1; y++)
        {
            uint32_t* row = mask + y * rs;

            row[sx0] |= mask0;

            for (int sx = sx0 + 1; sx < sx1; sx++)
                row[sx] = ~0;

            row[sx1] |= mask1;
        }
    }
}

void DOL::CreateBitMaskFromBlock(int w, int h, int sides, uint32_t mask[])
{
    int sw = MaskSize(w);
    int maskSize = sw * h;

    memset(mask, 0, sizeof(mask[0]) * maskSize);
    
    int x0 = w / 4;
    int x1 = x0 + w / 2;
    int y0 = h / 3;
    int y1 = y0 + h / 3;
    
    if (sides > 0)
        FillMaskBlock(sw, x0, x1, y0, y0, mask);
    if (sides > 1)
        FillMaskBlock(sw, x0, x0, y0, y1, mask);
    if (sides > 2)
        FillMaskBlock(sw, x0, x1, y1, y1, mask);
    if (sides > 3)
        FillMaskBlock(sw, x1, x1, y0, y1, mask);
}



// dirW is a w x h array of the fractional solid angle
// obscured in one of the four quadrants we're sweeping over.
void DOL::InitDirWFromBitMask(int w, int h, const uint32_t mask[], float dirW[])
{
    // mask is expected to have rows padded to whole numbers of uint32_ts
    int s = w * h;
    for (int i = 0; i < s * 4; i++)
        dirW[i] = 0.0f;
    
    int maskStride = MaskSize(w);
    
    for (int y = 0; y < h; y++)
        for (int j = 0; j < maskStride; j++)
        {
            uint32_t m = *mask++;
            
            if (m == 0)
                continue;
            
            float* block = dirW + y * w + 32 * j;
            
            for (int i = 0; i < 32; i++)
            {
                if (m & 1)
                {
                    block[i + 0 * s] = 1.0f;
                    block[i + 1 * s] = 1.0f;
                    block[i + 2 * s] = 1.0f;
                    block[i + 3 * s] = 1.0f;
                }
                
                m >>= 1;
            }
        }
}

namespace
{
    inline float QuadrantCombine(float a, float b)
    {
        // Find new quadrant estimate given same quad estimates from neighbour cells in
        // quad direction.
        // Simplest is to assume a and b cover half the quadrant each.
        // Need a=b=1 -> 1 to avoid light leaking.
        
    #if 0
        return 0.5f * (a + b) - 0.1f * a * b;
    #elif 0
        return 0.5f * (a + b);
    #else
        // if either leads to complete coverage, a + b - ab
        // if max k coverage from either, ka + kb - ksab? but would like a = b = 1 -> 1, get 2k - s = 1, s = 2k-1
        float k = 0.475f;
        float t = 0.975f;

        return t * (k * (a + b) - (2 * k - 1) * a * b);
    #endif
    }
}


void DOL::OcclusionSweep(int w, int h, float dirW[])
{
    int sliceStride = w * h;
    int lastRow = w * (h - 1);
    
    int sdx[4] = { +1, -1, +1, -1 };
    int sdy[4] = { +1, +1, -1, -1 };
    
    int sib[2] = { 0, w - 1 };
    int sie[2] = { w,    -1 };
    
    int rowStart[2] = { 0, lastRow };
    
    for (int sweep = 0; sweep < 4; sweep++)
    {
        const int dx = sdx[sweep];
        const int dy = sdy[sweep];
        
        const int sweepX = (sweep >> 0) & 1;
        const int sweepY = (sweep >> 1) & 1;

        const int ib = sib[sweepX];
        const int ie = sie[sweepX];

        float* row0 = (dirW + sweep * sliceStride) + rowStart[sweepY];
        
        // do first row -- only depends on itself
        for (int i = ib + dx; i != ie; i += dx)
            row0[i] += (1.0f - row0[i]) * QuadrantCombine(row0[i - dx], 0.0f);
            
        for (int j = 1; j < h; j++)
        {
            float* row1 = row0 + w * dy;
            
            // do first cell, only depends on prev row
            row1[ib] += (1.0f - row1[ib]) * QuadrantCombine(row0[ib], 0.0f);

            // rest of row -- prev cell and prev row
            for (int i = ib + dx; i != ie; i += dx)
                row1[i] += (1.0f - row1[i]) * QuadrantCombine(row1[i - dx], row0[i]);
            
            row0 = row1;
        }
    }
}


// --- 3D Volumes --------------------------------------------------------------

void DOL::CreateBitMask(int w, int h, int d, uint32_t seed, uint32_t mask[])
{
    int scale = (seed >> 4);
    seed &= 15;
    
    int count = (w + h + d) * (1 + scale);
    int ws = MaskSize(w);
    
    for (int i = 0; i < count; i++)
    {
        seed = NextSeed(seed);
        int r = seed & ~0x80000000;

        int s = r / w;
        int t = s / h;
        int u = t / d;
        
        int x = r - w * s;
        int y = s - h * t;
        int z = t - d * u;
        
        assert(x < w);
        assert(y < h);
        assert(z < d);
        
        int fx = x & 31;
        x >>= 5;
        
        int index = z * h * ws + y * ws + x;
        mask[index] |= 1 << fx;
    }
}

void DOL::InitDeltaFromBitMask(int w, int h, int d, const uint32_t mask[], cCellDelta3 delta[])
{
    cCellDelta3 maxDelta = { kMaxDelta3, kMaxDelta3, kMaxDelta3 };
    cCellDelta3 minDelta = { 0, 0, 0 };
    
    int s = w * h * d;
    for (int i = 0; i < s; i++)
        delta[i] = maxDelta;
    
    const int maskStride = MaskSize(w);
    
    for (int z = 0; z < d; z++)
        for (int y = 0; y < h; y++)
            for (int j = 0; j < maskStride * 32; j += 32)
            {
                uint32_t m = (*mask++);
                
                if (m == 0)
                    continue;
                
                cCellDelta3* block = delta + z * h * w + y * w + j;
                
                for (int i = 0; i < 32; i++)
                {
                    if (m & 1)
                        block[i] = minDelta;
                    
                    m >>= 1;
                }
            }
}

namespace
{
    inline bool Closer(cCellDelta3 d0, cCellDelta3 d1)
    {
        int dw0 = d0.x * d0.x + d0.y * d0.y + d0.z * d0.z;
        int dw1 = d1.x * d1.x + d1.y * d1.y + d1.z * d1.z;
        
        return dw0 > dw1;
    }

    void Danielsson2(int w, int h, cCellDelta3 delta[])
    {
        cCellDelta3* row0 = delta;
        cCellDelta3 cell;
            
        for (int j = 1; j < h; j++)
        {
            cCellDelta3* row1 = row0 + w;
            
            for (int i = 0; i < w; i++)
            {
                cell.x = row0[i].x;
                cell.y = row0[i].y - 1;
                cell.z = row0[i].z;
                
                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }
            
            for (int i = 1; i < w; i++)
            {
                cell.x = row1[i - 1].x - 1;
                cell.y = row1[i - 1].y;
                cell.z = row1[i - 1].z;
                
                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }

            for (int i = w - 2; i >= 0; i--)
            {

                cell.x = row1[i + 1].x + 1;
                cell.y = row1[i + 1].y;
                cell.z = row1[i + 1].z;
                
                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }
            
            row0 = row1;
        }
        
        row0 = delta + h * w - w;
        
        for (int j = h - 1; j > 0; j--)
        {
            cCellDelta3* row1 = row0 - w;
            
            for (int i = 0; i < w; i++)
            {
                cell.x = row0[i].x;
                cell.y = row0[i].y + 1;
                cell.z = row0[i].z;
                
                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }
            
            for (int i = 1; i < w; i++)
            {
                cell.x = row1[i - 1].x - 1;
                cell.y = row1[i - 1].y;
                cell.z = row1[i - 1].z;
                
                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }
            
            for (int i = w - 2; i >= 0; i--)
            {
                cell.x = row1[i + 1].x + 1;
                cell.y = row1[i + 1].y;
                cell.z = row1[i + 1].z;
                
                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }
            
            row0 = row1;
        }
    }
}

void DOL::Danielsson(int w, int h, int d, cCellDelta3 delta[])
{
    int sliceStride = w * h;
    
    cCellDelta3* slice0 = delta;
    cCellDelta3 cell;
    
    for (int k = 1; k < d; k++)
    {
        cCellDelta3* slice1 = slice0 + sliceStride;
        
        for (int i = 0; i < sliceStride; i++)
        {
            cell.x = slice0[i].x;
            cell.y = slice0[i].y;
            cell.z = slice0[i].z - 1;
            
            if (Closer(slice1[i], cell))
                slice1[i] = cell;
        }

        Danielsson2(w, h, slice1);
        
        slice0 = slice1;
    }

    slice0 = delta + sliceStride * (d - 1);
    
    for (int k = d - 2; k >= 0; k--)
    {
        cCellDelta3* slice1 = slice0 - sliceStride;
        
        for (int i = 0; i < sliceStride; i++)
        {
            cell.x = slice0[i].x;
            cell.y = slice0[i].y;
            cell.z = slice0[i].z + 1;
            
            if (Closer(slice1[i], cell))
                slice1[i] = cell;
        }
        
        Danielsson2(w, h, slice1);
        
        slice0 = slice1;
    }
}


void DOL::FastSweep(int w, int h, int d, cCellDelta3 delta[])
{
    // Fast sweep approach -- sweep each diagonal in turn.
    // C B
    // A x
    int sliceStride = w * h;
    int lastRow = w * (h - 1);
    
    int sdx[8] = { +1, -1, +1, -1, +1, -1, +1, -1 };
    int sdy[8] = { +1, +1, -1, -1, +1, +1, -1, -1 };
    int sdz[8] = { +1, +1, +1, +1, -1, -1, -1, -1 };
    
    int ib[2] = { 0, w - 1 };
    int ie[2] = { w,    -1 };
    
    int rowStart[2] = { 0, lastRow };
    
    for (int sweep = 0; sweep < 8; sweep++)
    {
        const int dx = sdx[sweep];
        const int dy = sdy[sweep];
        const int dz = sdz[sweep];
        
        const int sweepX = (sweep >> 0) & 1;
        const int sweepY = (sweep >> 1) & 1;
        const int sweepZ = (sweep >> 2) & 1;
        
        cCellDelta3* slice0 = delta + sweepZ * (d - 1) * sliceStride;
        cCellDelta3 cell;

        for (int k = 1; k < d; k++)
        {
            cCellDelta3* slice1 = slice0 + sliceStride * dz;
            
            cCellDelta3* row00 = slice0 + rowStart[sweepY];
            cCellDelta3* row01 = slice1 + rowStart[sweepY];
            
            for (int j = 1; j < h; j++)
            {
                cCellDelta3* row10 = row00 + w * dy;
                cCellDelta3* row11 = row01 + w * dy;
            
                for (int i = ib[sweepX]; i != ie[sweepX]; i += dx)
                {
                    cell.x = row10[i].x;
                    cell.y = row10[i].y;
                    cell.z = row10[i].z - dz;
                    
                    if (Closer(row11[i], cell))
                        row11[i] = cell;

                    cell.x = row01[i].x;
                    cell.y = row01[i].y - dy;
                    cell.z = row01[i].z;
                    
                    if (Closer(row11[i], cell))
                        row11[i] = cell;
                }

                for (int i = ib[sweepX] + dx; i != ie[sweepX]; i += dx)
                {
                    cell.x = row11[i - dx].x - dx;
                    cell.y = row11[i - dx].y;
                    cell.z = row11[i - dx].z;
                    
                    if (Closer(row11[i], cell))
                        row11[i] = cell;
                }
                
                row00 = row10;
                row01 = row11;
            }
            
            slice0 = slice1;
        }
    }
}



// Brute force -- for checking only!
void DOL::FindOptimalDeltas(int w, int h, int d, cCellDelta3 delta[])
{
    for (int k = 0; k < d; k++)
        for (int j = 0; j < h; j++)
            for (int i = 0; i < w; i++)
            {
                cCellDelta3* cell = delta + k * w * h + j * w + i;
                
                if (cell->x != 0 || cell->y != 0 || cell->z != 0)
                    continue;

                for (int z = 0; z < d; z++)
                    for (int y = 0; y < h; y++)
                        for (int x = 0; x < w; x++)
                        {
                            cCellDelta3* cell = delta + z * w * h + y * w + x;
                            
                            int d0 = (i - x) * (i - x) + (j - y) * (j - y) + (k - z) * (k - z); // d2 to i, j
                            int d1 = cell->x * cell->x + cell->y * cell->y + cell->z * cell->z; // current
                            
                            if (d0 < d1)
                            {
                                cell->x = i - x;
                                cell->y = j - y;
                                cell->z = k - z;
                            }
                        }                
            }
}

namespace
{
    void FillBlock(int rs, int ss, int x0, int x1, int y0, int y1, int z0, int z1, float* dirW[8])
    {
        for (int z = z0; z <= z1; z++)
            for (int y = y0; y <= y1; y++)
                for (int x = x0; x <= x1; x++)
                {
                    int i = z * ss + y * rs + x;
                    
                    for (int k = 0; k < 8; k++)
                        dirW[k][i] = 1.0f;
                }
    }

    void FillMaskBlock(int rs, int ss, int x0, int x1, int y0, int y1, int z0, int z1, uint32_t mask[])
    {
        for (int z = z0; z <= z1; z++)
        {
            uint32_t* sliceMask = mask + z * ss;

            FillMaskBlock(rs, x0, x1, y0, y1, sliceMask);
        }
    }
    
    void FillMaskBlockX(int rs, int ss, int x0, int x1, int y0, int y1, int z0, int z1, uint32_t mask[])
    {
        // e.g., 20 to 70:
        int sx0 = (x0 +  0) / 32;   // 0
        int sx1 = (x1 +  0) / 32;   // 2

        uint32_t mask0 = LeftMask (x0 - sx0 * 32);
        uint32_t mask1 = RightMask(x1 - sx1 * 32);

        if (sx0 == sx1)
        {
            mask0 = mask0 & mask1;
            mask1 = 0;
        }

        for (int z = z0; z <= z1; z++)
            for (int y = y0; y <= y1; y++)
            {
                uint32_t* row = mask + z * ss + y * rs;

                row[sx0] |= mask0;

                for (int sx = sx0 + 1; sx < sx1; sx++)
                    row[sx] = ~0;

                row[sx1] |= mask1;
            }
    }
}

void DOL::CreateBitMaskFromBlock(int w, int h, int d, int variant, uint32_t mask[])
{
    int sw = MaskSize(w);
    int rowStride   = sw;
    int sliceStride = sw * h;
    int volStride   = sw * h * d;

    memset(mask, 0, sizeof(mask[0]) * volStride);
    
    int x0 = w / 4;
    int x1 = x0 + w / 2;
    int y0 = h / 3;
    int y1 = y0 + h / 3;
    int z0 = d / 3;
    int z1 = z0 + d / 3;
    
    if (variant > 0) // -Z
        FillMaskBlock(rowStride, sliceStride, x0, x1, y0, y1, z0, z0, mask);
    if (variant > 1) // -X
        FillMaskBlock(rowStride, sliceStride, x0, x0, y0, y1, z0, z1, mask);
    if (variant > 2) // -Y
        FillMaskBlock(rowStride, sliceStride, x0, x1, y0, y0, z0, z1, mask);
    if (variant > 3) // +X
        FillMaskBlock(rowStride, sliceStride, x1, x1, y0, y1, z0, z1, mask);
    if (variant > 4) // +Y
        FillMaskBlock(rowStride, sliceStride, x0, x1, y1, y1, z0, z1, mask);
    if (variant > 5) // +Z
        FillMaskBlock(rowStride, sliceStride, x0, x1, y0, y1, z1, z1, mask);
}

void DOL::InitDirWFromBitMask(int w, int h, int d, const uint32_t mask[], float dirW[])
/// Initialise octants from bit mask.
{
    int s = w * h * d;
    for (int i = 0; i < s * 8; i++)
        dirW[i] = 0.0f;
    
    const int maskStride = MaskSize(w);
    
    for (int z = 0; z < d; z++)
        for (int y = 0; y < h; y++)
            for (int j = 0; j < maskStride * 32; j += 32)
            {
                uint32_t m = (*mask++);
                
                if (m == 0)
                    continue;
                
                float* block = dirW + z * h * w + y * w + j;
                
                for (int i = 0; i < 32; i++)
                {
                    if (m & 1)
                    {
                        for (int io = 0; io < 8; io++)
                            block[i + s * io] = 1.0f;
                    }
                    
                    m >>= 1;
                }
            }
}


namespace
{
    inline float OctantCombine(float a, float b, float c)
    {
    #if 1
        return 0.33333333f * (a + b + c - a * b * c);
    #elif 1
        return 0.33333333f * (a + b + c);
    #else

        // a + b + c - (ab + bc + ca) + abc   -> abc =s -> 3s - 3s^2 + s^3
        // if max k coverage from either, ka + kb + kc - (ab + bc + ac - abc) but would like a = b = 1 -> 1, get 2k - s = 1, s = 2k-1
        float k = 0.4f;
        float t = 0.975f;
        return t * (k * (a + b) - (2 * k - 1) * a * b);
    #endif
    }

    // ix = input occlusion, ox = output occlusion
    inline float OctantCombine(float oa, float ob, float oc, float ia, float ib, float ic)
    {
        float a = ia + (1.0f - ia) * oa;
        float b = ib + (1.0f - ib) * ob;
        float c = ic + (1.0f - ic) * oc;

        return 0.3333333f * (a + b + c - a * b * c);
    }

    void OcclusionSweep
    (
        int w,  int h,  int d,      // volume size
        int sx, int sy, int sz,     // signs of sweep direction

        float dirW[]                // in/out: seed occlusion, will be updated with swept occlusion
    )
    {
        const int dx = 1 - 2 * sx;
        const int dy = 1 - 2 * sy;
        const int dz = 1 - 2 * sz;
            
        const int ib = (w - 1) * sx;
        const int ie = w - sx * (w + 1);

        const int wh = w * h;   // slice stride
        const int whd = w * h * d;    // volume stride

        const int rowStart   = sy * (wh - w);
        const int sliceStart = sz * (whd - wh);

        float* slice0 = dirW + sliceStart;

        // do first slice -- only depends on itself
        {
            float* row0 = slice0 + rowStart;
            
            // as does its first row
            for (int i = ib + dx; i != ie; i += dx)
                row0[i] += (1.0f - row0[i]) * OctantCombine(row0[i - dx], 0.0f, 0.0f);
            
            for (int j = 1; j < h; j++)
            {
                float* row1 = row0 + w * dy;
                
                // do first cell, only depends on prev row
                row1[ib] += (1.0f - row1[ib]) * OctantCombine(0.0f, row0[ib], 0.0f);
                
                // rest of row -- prev cell and prev row
                for (int i = ib + dx; i != ie; i += dx)
                    row1[i] += (1.0f - row1[i]) * OctantCombine(row1[i - dx], row0[i], 0.0f);
                
                row0 = row1;
            }
        }
        
        // do remaining slices
        for (int k = 1; k < d; k++)
        {
            float* slice1 = slice0 + wh * dz;
            
            float* row00 = slice0 + rowStart;
            float* row01 = slice1 + rowStart;
            
            // do first cell of first row, only depends on prev slice
            row01[ib] += (1.0f - row01[ib]) * OctantCombine(0.0f, 0.0f, row00[ib]);
            
            // do first row -- only depends on itself and prev slice
            for (int i = ib + dx; i != ie; i += dx)
                row01[i] += (1.0f - row01[i]) * OctantCombine(row01[i - dx], 0.0f, row00[i]);

            for (int j = 1; j < h; j++)
            {
                float* row10 = row00 + w * dy;
                float* row11 = row01 + w * dy;
                
                // do first cell, only depends on prev row and prev slice
                row11[ib] += (1.0f - row11[ib]) * OctantCombine(0.0f, row01[ib], row10[ib]);
                
                // rest of row -- prev cell and prev row and prev slice
                for (int i = ib + dx; i != ie; i += dx)
                    row11[i] += (1.0f - row11[i]) * OctantCombine(row11[i - dx], row01[i], row10[i]);
                
                row00 = row10;
                row01 = row11;
            }
            
            slice0 = slice1;
        }
    }

    void OcclusionSweep
    (
        int w,  int h,  int d,      // volume size
        int sx, int sy, int sz,     // signs of sweep direction

        const float* inDirW,        // occlusion in this direction
        float*       dirW           // derived occlusion, not including original occlusion
    )
    {
        const int dx = 1 - 2 * sx;
        const int dy = 1 - 2 * sy;
        const int dz = 1 - 2 * sz;
            
        const int ib = (w - 1) * sx;
        const int ie = w - sx * (w + 1);

        const int wh = w * h;
        const int whd = w * h * d;

        const int rowStart   = sy * (wh - w);
        const int sliceStart = sz * (whd - wh);

        const float* islice0 = inDirW + sliceStart;
        float*       oslice0 = dirW + sliceStart;

        // do first slice -- only depends on itself
        {
            const float* irow0 = islice0 + rowStart;
            float*       orow0 = oslice0 + rowStart;
            
            // as does its first row
            for (int i = ib + dx; i != ie; i += dx)
                orow0[i] = OctantCombine(orow0[i - dx], 0.0f, 0.0f, irow0[i - dx], 0.0f, 0.0f);
            
            for (int j = 1; j < h; j++)
            {
                const float* irow1 = irow0 + w * dy;
                float*       orow1 = orow0 + w * dy;
                
                // do first cell, only depends on prev row
                orow1[ib] = OctantCombine(0.0f, orow0[ib], 0.0f, 0.0f, irow0[ib], 0.0f);
                
                // rest of row -- prev cell and prev row
                for (int i = ib + dx; i != ie; i += dx)
                    orow1[i] = OctantCombine(orow1[i - dx], orow0[i], 0.0f, irow1[i - dx], irow0[i], 0.0f);
                
                irow0 = irow1;
                orow0 = orow1;
            }
        }
        
        // do remaining slices
        for (int k = 1; k < d; k++)
        {
            const float* islice1 = islice0 + wh * dz;
            float*       oslice1 = oslice0 + wh * dz;
            
            const float* irow00 = islice0 + rowStart;
            float*       orow00 = oslice0 + rowStart;
            const float* irow01 = islice1 + rowStart;
            float*       orow01 = oslice1 + rowStart;
            
            // do first cell of first row, only depends on prev slice
            orow01[ib] = OctantCombine(0.0f, 0.0f, orow00[ib], 0.0f, 0.0f, irow00[ib]);
            
            // do first row -- only depends on itself and prev slice
            for (int i = ib + dx; i != ie; i += dx)
                orow01[i] = OctantCombine(orow01[i - dx], 0.0f, orow00[i], irow01[i - dx], 0.0f, irow00[i]);

            for (int j = 1; j < h; j++)
            {
                const float* irow10 = irow00 + w * dy;
                float*       orow10 = orow00 + w * dy;
                const float* irow11 = irow01 + w * dy;
                float*       orow11 = orow01 + w * dy;
                
                // do first cell, only depends on prev row and prev slice
                orow11[ib] = OctantCombine(0.0f, orow01[ib], orow10[ib], 0.0f, irow01[ib], irow10[ib]);
                
                // rest of row -- prev cell and prev row and prev slice
                for (int i = ib + dx; i != ie; i += dx)
                    orow11[i] = OctantCombine(orow11[i - dx], orow01[i], orow10[i], irow11[i - dx], irow01[i], irow10[i]);
                
                irow00 = irow10;
                orow00 = orow10;
                irow01 = irow11;
                orow01 = orow11;
            }
            
            islice0 = islice1;
            oslice0 = oslice1;
        }
    }
}

void DOL::OcclusionSweep(int w, int h, int d, float* dirW)
{
    for (int sweep = 0; sweep < 8; sweep++)
    {
        const int sx = (sweep >> 0) & 1;
        const int sy = (sweep >> 1) & 1;
        const int sz = (sweep >> 2) & 1;

        ::OcclusionSweep(w, h, d, sx, sy, sz, dirW + sweep * w * h * d);
    }
}

void DOL::OcclusionSweep(int w, int h, int d, const float* occDirW, float* dirW)
{
    int whd = w * h * d;

    for (int sweep = 0; sweep < 8; sweep++)
    {
        const int sx = (sweep >> 0) & 1;
        const int sy = (sweep >> 1) & 1;
        const int sz = (sweep >> 2) & 1;

        ::OcclusionSweep(w, h, d, sx, sy, sz, occDirW + sweep * whd, dirW + sweep * whd);
    }
}
