//------------------------------------------------------------------------------
// Purpose: Distance Field and Occlusion Generation
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#include "DistanceOcclusionLib.h"

#include <math.h>
#include <string.h>

using namespace DOL;

#if defined(DOS_DEBUG) || defined(DEBUG)
    #include <assert.h>
#else
    #define assert(e) ((void) 0)
#endif

#if !defined(__cplusplus) || (__cplusplus <= 201100)
    #define constexpr const
#endif

// Declarations

namespace
{
    inline int RNG(uint32_t& seed, int limit)
    {
        seed = uint32_t(seed * uint64_t(1103515245)) + 12345;
        return int((seed * (uint64_t(limit))) >> 32);
    }
}

// Convenience macro for iterating over set bits (and also testing different strategies)
#ifdef _MSC_VER
    #include <intrin.h>
    uint32_t TrailingZeroes32(uint32_t mask)
    {
        unsigned long index;
        _BitScanForward(&index, mask);
        return index;
    }
#elif defined(__GNUC__)
    #define TrailingZeroes32 __builtin_ctz
#else
    // Hash table approach is actually comparable in speed to intrinsics, faster than hierarchical bit tricks
    constexpr int kTrailingZeroesHashTable[37] = { 32, 0, 1, 26, 2, 23, 27, 0, 3, 16, 24, 30, 28, 11, 0, 13, 4, 7, 17, 0, 25, 22, 31, 15, 29, 10, 12, 6, 0, 21, 14, 9, 5, 20, 8, 19, 18 };

    int TrailingZeroes32(uint32_t x)
    {
        return kTrailingZeroesHashTable[(x & -x) % 37];
    }
#endif

#define ITER_SET_BITS_BEGIN(MASK)           \
    while (MASK)                            \
    {                                       \
        int i = TrailingZeroes32(MASK);     \
        MASK &= MASK - 1;                   \
                                            \
        if (j * 32 + i >= w)                \
            break;                          \

#define ITER_SET_BITS_BEGIN_BF(MASK)        \
    for (int i = 0; i < 32; i++)            \
        if (MASK & (1 << i))                \
        {                                   \
            if (j * 32 + i >= w)            \
                break;                      \

#define ITER_SET_BITS_END }




// --- 2D Images ---------------------------------------------------------------

uint32_t* DOL::CreateBitMask(int w, int h)
{
    size_t maskSize = MaskSize(w, h);
    uint32_t* mask = new uint32_t[maskSize];

    memset(mask, 0, maskSize * sizeof(uint32_t));
    
    return mask;
}

void DOL::DestroyBitMask(uint32_t*& mask)
{
    delete[] mask;
    mask = 0;
}

void DOL::BitMaskAddPoints(int w, int h, int count, uint32_t seed, uint32_t mask[])
{
    int ws = MaskSize(w);
    
    for (int i = 0; i < count; i++)
    {
        int x = RNG(seed, w);
        int y = RNG(seed, h);

        assert(x < w);
        assert(y < h);
        
        int fx = x & 31;
        x >>= 5;
        
        int index = y * ws + x;
        mask[index] |= (1 << fx);
    }
}

namespace
{
    uint32_t LeftMask(int n)
    {
        assert(n < 32);
        uint32_t bit = (1 << n);
        return ~(bit - 1);
    }

    uint32_t RightMask(int n)
    {
        assert(n < 32);
        uint32_t bit = 1 << n;
        return bit + (bit - 1); // Inclusive
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

void DOL::BitMaskAddBlock(int w, int h, int sides, uint32_t mask[])
{
    int x0 = w / 4;
    int x1 = x0 + w / 2 - 1;
    int y0 = h / 3;
    int y1 = y0 + h / 3 - 1;

    int sw = MaskSize(w);

    if (sides > 4)
        return FillMaskBlock(sw, x0, x1, y0, y1, mask);

    if (sides > 0)
        FillMaskBlock(sw, x0, x1, y0, y0, mask);
    if (sides > 1)
        FillMaskBlock(sw, x0, x0, y0, y1, mask);
    if (sides > 2)
        FillMaskBlock(sw, x0, x1, y1, y1, mask);
    if (sides > 3)
        FillMaskBlock(sw, x1, x1, y0, y1, mask);
}

namespace
{
    template<class T> void InitDistancesFromBitMask(int w, int h, const uint32_t mask[], T distances[], T maxD)
    {
        // mask is expected to have rows padded to whole numbers of uint32_ts
        int s = w * h;
        for (int i = 0; i < s; i++)
            distances[i] = maxD;

        int maskStride = MaskSize(w);
        
        for (int y = 0; y < h; y++)
        {
            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m = *mask++;
                
                if (m == 0)
                    continue;
                
                T* block = distances + y * w + 32 * j;
                
                ITER_SET_BITS_BEGIN(m)
                    block[i] = 0;
                ITER_SET_BITS_END
            }
        }
    }
}

void DOL::InitDistancesFromBitMask(int w, int h, const uint32_t mask[], int16_t distances[], int16_t maxD)
{
    return ::InitDistancesFromBitMask<int16_t>(w, h, mask, distances, maxD);
}
void DOL::InitDistancesFromBitMask(int w, int h, const uint32_t mask[], int32_t distances[], int32_t maxD)
{
    return ::InitDistancesFromBitMask<int32_t>(w, h, mask, distances, maxD);
}
void DOL::InitDistancesFromBitMask(int w, int h, const uint32_t mask[], float distances[], float maxD)
{
    return ::InitDistancesFromBitMask<float>(w, h, mask, distances, maxD);
}

void DOL::InitDeltasFromBitMask(int w, int h, const uint32_t mask[], cCellDelta2 deltas[], bool invert)
{
    const cCellDelta2 maxDelta = { kMaxDelta2, kMaxDelta2 };
    const cCellDelta2 minDelta = { 0, 0 };
    
    for (int i = 0; i < w * h; i++)
        deltas[i] = maxDelta;
    
    const int      maskStride = MaskSize(w);
    const uint32_t maskInvert = invert ? ~uint32_t(0) : 0;

    for (int y = 0; y < h; y++)
        for (int j = 0; j < maskStride * 32; j += 32)
        {
            uint32_t m = (*mask++) ^ maskInvert;
            
            if (m == 0)
                continue;
            
            cCellDelta2* block = deltas + y * w + j;
            const int n = (w - j) >= 32 ? 32 : (w - j);

            for (int i = 0; i < n; i++) // XXX
            {
                if (m & 1)
                    block[i] = minDelta;
                
                m >>= 1;
            }
        }
}

namespace
{
    // Utilities for setting deltas for cells adjacent to the interior/exterior
    // boundary, in units of half a cell width. E.g., a cell immediately to the
    // left of the boundary will have a delta of (+1, 0), reflecting that the
    // boundary is half a cell to right of the cell's centre.
    // By initialising the SDF this way, we can calculate both interior and
    // exterior distances in a single pass, and distances are consistent across
    // the border. (The standard techique of calculating the exterior distances,
    // then inverting the problem to calculate interior distances, and combining
    // the two results, is both 2x as expensive, and leads to inconsistent
    // distances on border cells. The border effectively has a width of -ve 1.
    
    constexpr cCellDelta2 kBH0 = { 0, +1 };  // level set is half a cell up
    constexpr cCellDelta2 kBH1 = { 0, -1 };  // level set is half a cell down
    constexpr cCellDelta2 kBW0 = { +1, 0 };  // level set is half a cell right
    constexpr cCellDelta2 kBW1 = { -1, 0 };  // level set is half a cell left

#ifdef USE_BORDER_REFERENCE
    // These routines treat each border case separately, in turn, for clarity.
    // They are useful both for understanding the logic, but also for sanity-
    // checking the more complex all-in-one routine. (See SetBorderCellsCombo.)
    void SetBorderCellsH(int w, int h, const uint32_t* const mask, cCellDelta2 deltas[])
    {
        // Find horizontal edges...
        //   x        .
        //   .   or   x
        // and set the cells to either side accordingly

        const int maskStride = MaskSize(w);

        for (int y = 0; y < h - 1; y++)
        {
            const uint32_t* maskRow0   = mask   + (y + 0) * maskStride;
            const uint32_t* maskRow1   = mask   + (y + 1) * maskStride;
            cCellDelta2*    deltasRow0 = deltas + (y + 0) * w;
            cCellDelta2*    deltasRow1 = deltas + (y + 1) * w;

            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m0 = maskRow0[j];
                uint32_t m1 = maskRow1[j];
                uint32_t mx = m0 ^ m1;  // 1 on change

                if (mx)
                {
                    cCellDelta2* block0 = deltasRow0 + j * 32;
                    cCellDelta2* block1 = deltasRow1 + j * 32;

                    for (int i = 0; i < 32; i++)
                        if ((mx & (1 << i)))
                        {
                            if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                                break;

                            assert(j * 32 + i >= 0);
                            assert(j * 32 + i < w);

                            block0[i] = kBH0;
                            block1[i] = kBH1;
                        }
                }
            }
        }
    }

    void SetBorderCellsW(int w, int h, const uint32_t* const mask, cCellDelta2 deltas[])
    {
        // Find vertical edges...
        //   .|x   or   x|.
        // and set the cells to either side accordingly

        const int maskStride = MaskSize(w);

        for (int y = 0; y < h; y++)
        {
            const uint32_t* maskRow   = mask   + y * maskStride;
            cCellDelta2*    deltasRow = deltas + y * w;

            uint32_t lastBit = maskRow[0] & 1; // repeat border cell, avoids extra logic to loop from 1 rather than 0.

            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m0 = maskRow[j];
                uint32_t m1 = (m0 << 1) | lastBit;
                lastBit = m0 >> 31;

                // Effectively bit 'i' indexes cell i (m0) and cell i-1 (m1). It's
                // done this way to avoid having to look ahead one word up front.

                uint32_t mx = m0 ^ m1;

                if (mx)
                {
                    cCellDelta2* block = deltasRow + j * 32;

                    for (int i = 0; i < 32; i++)
                        if ((mx & (1 << i)))
                        {
                            if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                                break;

                            assert(j * 32 + i > 0);
                            assert(j * 32 + i < w);

                            block[i - 1] = kBW0;
                            block[i + 0] = kBW1;
                        }
                }
            }
        }
    }

    constexpr cCellDelta2 kBC00 = { -1, +1 };
    constexpr cCellDelta2 kBC01 = { +1, -1 };
    constexpr cCellDelta2 kBC10 = { +1, +1 };
    constexpr cCellDelta2 kBC11 = { -1, -1 };

    void SetBorderCellsC(int w, int h, const uint32_t* const mask, cCellDelta2 deltas[])
    {
        // Look for diagonal borders affecting corners, i.e.,
        //  . x       x .
        //  x .   or  . x
        // and set the diagonal cell deltas to +-1, +-1 as appropriate.
        // These may be overwritten by horizontal/vertical deltas later, as they
        // will always be closer.

        const int maskStride = MaskSize(w);

        for (int y = 0; y < h - 1; y++)
        {
            const uint32_t* maskRow0   = mask   + (y + 0) * maskStride;
            const uint32_t* maskRow1   = mask   + (y + 1) * maskStride;
            cCellDelta2*    deltasRow0 = deltas + (y + 0) * w;
            cCellDelta2*    deltasRow1 = deltas + (y + 1) * w;

            uint32_t lastBit0 = maskRow1[0] & 1;    // pick to avoid detection on i=j=0 without explicit check
            uint32_t lastBit1 = maskRow0[0] & 1;

            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m00 = maskRow0[j];
                uint32_t m01 = maskRow1[j];

                uint32_t m10 = (m00 << 1) | lastBit0;
                uint32_t m11 = (m01 << 1) | lastBit1;

                lastBit0 = m00 >> 31;
                lastBit1 = m01 >> 31;

                uint32_t mc0 = m00 ^ m11;   // forward diagonal
                uint32_t mc1 = m01 ^ m10;   // backward diagonal

                if (mc0 | mc1)
                {
                    cCellDelta2* block[2] = { deltasRow0 + j * 32, deltasRow1 + j * 32 };

                    for (int i = 0; i < 32; i++)
                    {
                        if ((mc0 & (1 << i)))
                        {
                            if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                                break;

                            assert(j * 32 + i > 0);
                            assert(j * 32 + i < w);

                            block[0][i + 0] = kBC00;
                            block[1][i - 1] = kBC01;
                        }
                        if ((mc1 & (1 << i)) && (j * 32 + i < w))
                        {
                            if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                                break;

                            assert(j * 32 + i > 0);
                            assert(j * 32 + i < w);

                            block[0][i - 1] = kBC10;
                            block[1][i + 0] = kBC11;
                        }
                    }
                }
            }
        }
    }

    void SetBorderCellsC2(int w, int h, const uint32_t* const mask, cCellDelta2 deltas[])
    {
        // Smarter version that figures out which particular cell
        // of the diagonal should be set. Reduces writes by 2x.

        const int maskStride = MaskSize(w);

        for (int y = 0; y < h - 1; y++)
        {
            const uint32_t* maskRow0   = mask   + (y + 0) * maskStride;
            const uint32_t* maskRow1   = mask   + (y + 1) * maskStride;
            cCellDelta2*    deltasRow0 = deltas + (y + 0) * w;
            cCellDelta2*    deltasRow1 = deltas + (y + 1) * w;

            uint32_t lastBit0 = maskRow0[0] & 1;
            uint32_t lastBit1 = maskRow1[0] & 1;

            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m00 = maskRow0[j];
                uint32_t m01 = maskRow1[j];

                uint32_t m10 = (m00 << 1) | lastBit0;
                uint32_t m11 = (m01 << 1) | lastBit1;

                lastBit0 = m00 >> 31;
                lastBit1 = m01 >> 31;

                // we have a corner cell if both verticals (or horizontals) feature one flip
                uint32_t mv0 = m00 ^ m01;  // 1 on vert change
                uint32_t mv1 = m10 ^ m11;  // 1 on vert change
                uint32_t mh0 = m00 ^ m10;  // 1 on change

                uint32_t mc = (mv0 ^ mv1);

                if (mc)
                {
                    cCellDelta2* block[2] = { deltasRow0 + j * 32, deltasRow1 + j * 32 };

                    for (int i = 0; i < 32; i++)
                        if (mc & (1 << i))
                        {
                            if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                                break;

                            assert(j * 32 + i > 0);
                            assert(j * 32 + i < w);

                            int dx = (mv1 >> i) & 1;
                            int dy = (mh0 >> i) & 1;

                            block[dy][i + dx - 1].x = 1 - 2 * dx;
                            block[dy][i + dx - 1].y = 1 - 2 * dy;
                        }
                }
            }
        }
    }

#else

    void SetBorderCellsCombo(int w, int h, const uint32_t* const mask, cCellDelta2 deltas[])
    {
        // Rolls everything into one pass
        const int maskStride = MaskSize(w);

        for (int y = 0; y < h - 1; y++)
        {
            const uint32_t* maskRow0   = mask   + (y + 0) * maskStride;
            const uint32_t* maskRow1   = mask   + (y + 1) * maskStride;
            cCellDelta2*    deltasRow0 = deltas + (y + 0) * w;
            cCellDelta2*    deltasRow1 = deltas + (y + 1) * w;

            uint32_t lastBit0 = maskRow0[0] & 1;    // pick to avoid vertical i=0 check
            uint32_t lastBit1 = maskRow1[0] & 1;

            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m00 = maskRow0[j];
                uint32_t m01 = maskRow1[j];
                uint32_t m10 = (m00 << 1) | lastBit0;
                uint32_t m11 = (m01 << 1) | lastBit1;

                lastBit0 = m00 >> 31;
                lastBit1 = m01 >> 31;

                cCellDelta2* block[2] = { deltasRow0 + j * 32, deltasRow1 + j * 32 };

                // Horizontal
                uint32_t ey = m00 ^ m01;  // 1 on change

                // Corner. These are much rarer, and we must check to ensure we don't
                // overwrite edge deltas, as they are larger.
                uint32_t ex = m00 ^ m10;
                uint32_t eyx = m10 ^ m11;
                uint32_t bc = ey ^ eyx;

                // The order here is important as iterating over the masks can destroy them

                // if (ex | ey | bc) early out currently not worth it
                ITER_SET_BITS_BEGIN(bc)
                    assert(j * 32 + i > 0);

                    int dx = (eyx >> i) & 1;
                    int dy = (ex  >> i) & 1;

                    cCellDelta2& corner = block[dy][i + dx - 1];

                    if (corner.x != kMaxDelta2)
                        continue;

                    corner.x = 1 - 2 * dx;
                    corner.y = 1 - 2 * dy;
                ITER_SET_BITS_END

                ITER_SET_BITS_BEGIN(ex)
                    assert(j * 32 + i > 0);
                    block[0][i - 1] = kBW0;
                    block[0][i + 0] = kBW1;
                ITER_SET_BITS_END

                ITER_SET_BITS_BEGIN(ey)
                    block[0][i] = kBH0;
                    block[1][i] = kBH1;
                ITER_SET_BITS_END
            }
        }

        // Finish last vertical row
        {
            const int y = h - 1;
            const uint32_t* maskRow   = mask   + y * maskStride;
            cCellDelta2*    deltasRow = deltas + y * w;

            uint32_t lastBit = maskRow[0] & 1; // repeat border cell, avoids extra logic to loop from 1 rather than 0.

            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m0 = maskRow[j];
                uint32_t m1 = (m0 << 1) | lastBit;
                lastBit = m0 >> 31;

                uint32_t mvx = m0 ^ m1;

                if (mvx)
                {
                    cCellDelta2* block = deltasRow + j * 32;

                    ITER_SET_BITS_BEGIN(mvx)
                        assert(j * 32 + i > 0);
                        block[i - 1] = kBW0;
                        block[i + 0] = kBW1;
                    ITER_SET_BITS_END
                }
            }
        }
    }
#endif
}


void DOL::InitBorderDeltasFromBitMask(int w, int h, const uint32_t* const mask, cCellDelta2 deltas[])
{
    const cCellDelta2 maxDelta = { kMaxDelta2, kMaxDelta2 };

    int s = w * h;
    for (int i = 0; i < s; i++)
        deltas[i] = maxDelta;

#ifdef USE_BORDER_REFERENCE
    SetBorderCellsC(w, h, mask, deltas);
    SetBorderCellsH(w, h, mask, deltas);
    SetBorderCellsW(w, h, mask, deltas);
#else
    SetBorderCellsCombo(w, h, mask, deltas);
#endif
}

namespace
{
    /// Returns true if d1 is closer than d0
    inline bool Closer(cCellDelta2 d0, cCellDelta2 d1)
    {
        int dw0 = d0.x * d0.x + d0.y * d0.y;
        int dw1 = d1.x * d1.x + d1.y * d1.y;
        
        return dw0 > dw1;
    }
}

void DOL::FastSweep(int w, int h, cCellDelta2 deltas[], int cw)
{
    // Fast sweep approach -- sweep each diagonal in turn.
    // This is really the canonical sweeping method -- each sweep calculates
    // distances over a quadrant, as for each cell in the sweep, you can 
    // guarantee all cells in the corresponding quadrant defined by that cell
    // and the start point have already been visited.
    // Other approaches such as Danielsson and Chamfer combine sweeps and vary
    // the inspected neighbours to trade speed vs. accuracy.
    // 
    // See, e.g., "A fast sweeping method for Eikonal equations" Zhao 2004.
    
    const int lastRow = w * (h - 1);
    cCellDelta2 cell;

    for (int sweep = 0; sweep < 4; sweep++)
    {
        const int sx = (sweep >> 0) & 1;
        const int sy = (sweep >> 1) & 1;

        const int dx = 1 - 2 * sx;
        const int dy = 1 - 2 * sy;
        const int cx = -cw * dx;
        const int cy = -cw * dy;

        const int ib =     sx * (w - 1);
        const int ie = w - sx * (w + 1);

        cCellDelta2* row0 = deltas + sy * lastRow;

        for (int j = 1; j < h; j++)
        {
            cCellDelta2* row1 = row0 + w * dy;

            for (int i = ib; i != ie; i += dx)
            {
                // Testing showed early out check for c = 0, 0 not worth it
                // even with large contiguous areas of interior cells. Actively
                // harmful in other cases due to branch overhead.
                cell.x = row0[i].x;
                cell.y = row0[i].y + cy;

                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }

            for (int i = ib + dx; i != ie; i += dx)
            {
                cell.x = row1[i - dx].x + cx;
                cell.y = row1[i - dx].y;

                if (Closer(row1[i], cell))
                    row1[i] = cell;

                cell.x = row0[i - dx].x + cx;
                cell.y = row0[i - dx].y + cy;

                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }

            row0 = row1;
        }
    }
}

// Danielsson is similar to the quadrant sweep approach (FastSweep), except it
// wraps alternating x scans into one vertical sweep.

void DOL::Danielsson(int w, int h, cCellDelta2 deltas[], int cw)
{
    // store dx^2, dy^2 to nearest point. distance = sqrt(dx^2 + dy^2).
    // Much more accurate, doesn't over-estimate like Chamfer.
    cCellDelta2* row0 = deltas;
    cCellDelta2 cell;

    for (int j = 1; j < h; j++)
    {
        cCellDelta2* row1 = row0 + w;
        
        for (int i = 0; i < w; i++)
        {
            cell.x = row0[i].x;
            cell.y = row0[i].y - cw;

            if (Closer(row1[i], cell))
                row1[i] = cell;
        }
        
        for (int i = 1; i < w; i++)
        {
            cell.x = row1[i - 1].x - cw;
            cell.y = row1[i - 1].y;

            if (Closer(row1[i], cell))
                row1[i] = cell;
        }

        for (int i = w - 2; i >= 0; i--)
        {
            cell.x = row1[i + 1].x + cw;
            cell.y = row1[i + 1].y;
            
            if (Closer(row1[i], cell))
                row1[i] = cell;
        }

        row0 = row1;
    }

    row0 = deltas + h * w - w;

    for (int j = h - 1; j > 0; j--)
    {
        cCellDelta2* row1 = row0 - w;
        
        for (int i = 0; i < w; i++)
        {
            cell.x = row0[i].x;
            cell.y = row0[i].y + cw;
            
            if (Closer(row1[i], cell))
                row1[i] = cell;
        }
        
        for (int i = 1; i < w; i++)
        {
            cell.x = row1[i - 1].x - cw;
            cell.y = row1[i - 1].y;
            
            if (Closer(row1[i], cell))
                row1[i] = cell;
        }

        for (int i = w - 2; i >= 0; i--)
        {
            cell.x = row1[i + 1].x + cw;
            cell.y = row1[i + 1].y;
            
            if (Closer(row1[i], cell))
                row1[i] = cell;
        }
        
        row0 = row1;
    }
}

namespace
{
    void JumpFlood(int step, int w, int h, cCellDelta2 deltasIn[], cCellDelta2 deltasOut[], int cw)
    {
        const cCellDelta2* row    = deltasIn;
        cCellDelta2*       rowOut = deltasOut;

        for (int j = 0; j < h; j++)
        {
            for (int i = 0; i < w; i++)
            {
                cCellDelta2 minCell = row[i]; 
                
                if (i - step >= 0)
                {
                    cCellDelta2 cell = row[i - step];
                    cell.x -= cw * step;
                    
                    if (Closer(minCell, cell))
                        minCell = cell;

                    if (j - step >= 0)
                    {
                        cCellDelta2 cell = row[i - step - step * w];
                        cell.x -= cw * step;
                        cell.y -= cw * step;

                        if (Closer(minCell, cell))
                            minCell = cell;
                    }

                    if (j + step < h)
                    {
                        cCellDelta2 cell = row[i - step + step * w];
                        cell.x -= cw * step;
                        cell.y += cw * step;

                        if (Closer(minCell, cell))
                            minCell = cell;
                    }
                }

                if (i + step < w)
                {
                    cCellDelta2 cell = row[i + step];
                    cell.x += cw * step;

                    if (Closer(minCell, cell))
                        minCell = cell;

                    if (j - step >= 0)
                    {
                        cCellDelta2 cell = row[i + step - step * w];
                        cell.x += cw * step;
                        cell.y -= cw * step;

                        if (Closer(minCell, cell))
                            minCell = cell;
                    }

                    if (j + step < h)
                    {
                        cCellDelta2 cell = row[i + step + step * w];
                        cell.x += cw * step;
                        cell.y += cw * step;

                        if (Closer(minCell, cell))
                            minCell = cell;
                    }
                }

                if (j - step >= 0)
                {
                    cCellDelta2 cell = row[i - step * w];
                    cell.y -= cw * step;

                    if (Closer(minCell, cell))
                        minCell = cell;
                }

                if (j + step < h)
                {
                    cCellDelta2 cell = row[i + step * w];
                    cell.y += cw * step;

                    if (Closer(minCell, cell))
                        minCell = cell;
                }

                rowOut[i] = minCell; 
            }
            
            row    += w;
            rowOut += w;
        }
    }

    int FloorPow2(int v)
    {
        v |= v >> 1;
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;

        return (v + 1) >> 1;
    }
}

void DOL::JumpFlood(int w, int h, cCellDelta2 deltas[], int cw)
{
    cCellDelta2* temp = new cCellDelta2[w * h];
    cCellDelta2* d0 = deltas;
    cCellDelta2* d1 = temp;

    int step = FloorPow2((w >= h ? w : h) - 1);

    for ( ; step > 0; step /= 2)
    {
        ::JumpFlood(step, w, h, d0, d1, cw);

        cCellDelta2* dt = d0;
        d0 = d1;
        d1 = dt;
    }

    if (d0 != deltas)
        ::JumpFlood(1, w, h, d0, deltas, cw);    // we could copy, but might as well do another round

    delete[] temp;
}


namespace
{
    // G. Borgefors.
    // Distance transformations in digital images.

    template<class T> int Chamfer(int w, int h, T distances[])
    {
        // according to Borgefors, using these values leads to more
        // stable and accurate results than 1, sqrt(2) with double precision.
        // See Table 1 in Borgefors '86. The resulting distances must be
        // divided by d0 once the algorithm is complete.

        const T d0 = 3;
        const T d1 = 4;

        T* row0 = distances;

        for (int x = 1; x < w; x++)             // restricted mask for y=0
            if (row0[x] > row0[x - 1] + d0)
                row0[x] = row0[x - 1] + d0;

        for (int y = 1; y < h; y++)
        {
            T* row1 = row0 + w;

            if (row1[0] > row0[0] + d0)         // restricted mask for x=0
                row1[0] = row0[0] + d0;
            if (row1[0] > row0[1] + d1)
                row1[0] = row0[1] + d1;

            for (int x = 1; x < w - 1; x++)     // full mask
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

            if (row1[w - 1] > row0[w - 2] + d1) // restricted mask for x=w-1
                row1[w - 1] = row0[w - 2] + d1;
            if (row1[w - 1] > row0[w - 1] + d0)
                row1[w - 1] = row0[w - 1] + d0;

            if (row1[w - 1] > row1[w - 2] + d0)
                row1[w - 1] = row1[w - 2] + d0;

            row0 = row1;
        }

        row0 = distances + h * w - w;

        for (int x = w - 2; x >= 0; x--)        // restricted reverse mask for y=h-1
            if (row0[x] > row0[x + 1] + d0)
                row0[x] = row0[x + 1] + d0;

        for (int y = h - 2; y >= 0; y--)
        {
            T* row1 = row0 - w;

            if (row1[w - 1] > row0[w - 1] + d0) // restricted mask for x=w-1
                row1[w - 1] = row0[w - 1] + d0;
            if (row1[w - 1] > row0[w - 2] + d1)
                row1[w - 1] = row0[w - 2] + d1;

            for (int x = w - 2; x > 0; x--)     // full reverse mask
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

            if (row1[0] > row0[1] + d1)         // restricted mask for x=0
                row1[0] = row0[1] + d1;
            if (row1[0] > row0[0] + d0)
                row1[0] = row0[0] + d0;

            if (row1[0] > row1[1] + d0)
                row1[0] = row1[1] + d0;

            row0 = row1;
        }
        
        return d0;
    }
}

int DOL::Chamfer(int w, int h, int16_t distances[])
{
    return ::Chamfer<int16_t>(w, h, distances);
}

int DOL::Chamfer(int w, int h, int32_t distances[])
{
    return ::Chamfer<int32_t>(w, h, distances);
}

namespace
{
    // "Distance Transforms of Sampled Functions", Felzenszwalb and Huttenlocher
    // Uses parabolas to track distance. Is still O(n), and allows float-valued
    // initial per-cell distance estimates, but the parabolas are still centred
    // on the cell, which limits accuracy for higher-than-cell-resolution
    // boundaries. Poor locality, is noticeably slower than the other O(n)
    // algorithms unless everything is in cache.
    void Felzenszwalb(float f[], int n, int v[/* n */], float z[/*n + 1*/])    // 1D version
    {
        int k = 0;

        v[0] = 0;
        z[0] = -FLT_MAX;
        z[1] = +FLT_MAX;

        for (int q = 1; q < n; q++)
        {
            float s  = ((f[q] + (q * q)) - (f[v[k]] + (v[k] * v[k]))) / (2 * q - 2 * v[k]);

            while (s <= z[k])
            {
                k--;
                s  = ((f[q] + (q * q)) - (f[v[k]] + (v[k] * v[k]))) / (2 * q - 2 * v[k]);
            }

            k++;
            v[k] = q;
            z[k] = s;
            z[k+1] = +FLT_MAX;
        }
        
        k = 0;
        for (int q = 0; q < n; q++)
        {
            while (z[k + 1] < q)
                k++;

            f[q] = (q - v[k]) * (q - v[k]) + f[v[k]];
        }
    }
}

void DOL::Felzenszwalb(int w, int h, float distances[])
{
    int whMax = w >= h ? w : h;
    float* f = new float[whMax];
    int*   v = new int[whMax];
    float* z = new float[whMax + 1];

    // Transform along columns
    for (int x = 0; x < w; x++)
    {
        for (int y = 0; y < h; y++)
            f[y] = distances[x + y * w];

        ::Felzenszwalb(f, h, v, z);

        for (int y = 0; y < h; y++)
            distances[x + y * w] = f[y];
    }
    
    // transform along rows
    for (int y = 0; y < h; y++)
        ::Felzenszwalb(distances + y * w, w, v, z);
    
    delete[] z;
    delete[] v;
    delete[] f;
}



// Brute force -- for checking only!
void DOL::BruteForce(int w, int h, cCellDelta2 deltas[], int cw)
{
    for (int i = 0; i < h; i++)
        for (int j = 0; j < w; j++)
        {
            cCellDelta2* cell = deltas + i * w + j;

            if (abs(cell->x) >= cw || abs(cell->y) >= cw)   // non-boundary cell?
                continue;

            // For each occupied cell, iterate over all other cells and check for minimal distance

            for (int y = 0; y < h; y++)
                for (int x = 0; x < w; x++)
                {
                    cCellDelta2 d2 = { int16_t((j - x) * cw + cell->x), int16_t((i - y) * cw + cell->y) };
                    cCellDelta2* cell2 = deltas + y * w + x;

                    if (Closer(*cell2, d2))
                        *cell2 = d2;
                }
        }
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
    const int sliceStride = w * h;
    const int lastRow = w * (h - 1);
    
    for (int sweep = 0; sweep < 4; sweep++)
    {
        const int sx = (sweep >> 0) & 1;
        const int sy = (sweep >> 1) & 1;

        const int dx = 1 - 2 * sx;
        const int dy = 1 - 2 * sy;

        const int ib =     sx * (w - 1);
        const int ie = w - sx * (w + 1);

        float* row0 = (dirW + sweep * sliceStride) + sy * lastRow;
        
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

uint32_t* DOL::CreateBitMask(int w, int h, int d)
{
    size_t maskSize = MaskSize(w, h, d);
    uint32_t* mask = new uint32_t[maskSize];

    memset(mask, 0, maskSize * sizeof(uint32_t));

    return mask;
}

void DOL::BitMaskAddPoints(int w, int h, int d, int count, uint32_t seed, uint32_t mask[])
{
    int ws = MaskSize(w);
    
    for (int i = 0; i < count; i++)
    {
        int x = RNG(seed, w);
        int y = RNG(seed, h);
        int z = RNG(seed, d);
        
        assert(x < w);
        assert(y < h);
        assert(z < d);
        
        int fx = x & 31;
        x >>= 5;
        
        int index = z * h * ws + y * ws + x;
        mask[index] |= 1 << fx;
    }
}

namespace
{
    void FillMaskBlock(int rs, int ss, int x0, int x1, int y0, int y1, int z0, int z1, uint32_t mask[])
    {
        for (int z = z0; z <= z1; z++)
        {
            uint32_t* sliceMask = mask + z * ss;

            FillMaskBlock(rs, x0, x1, y0, y1, sliceMask);
        }
    }
}

void DOL::BitMaskAddBlock(int w, int h, int d, int sides, uint32_t mask[])
{
    int sw = MaskSize(w);
    int rowStride   = sw;
    int sliceStride = sw * h;

    int x0 = w / 4;
    int x1 = x0 + w / 2 - 1;
    int y0 = h / 3;
    int y1 = y0 + h / 3 - 1;
    int z0 = d / 3;
    int z1 = z0 + d / 3 - 1;

    if (sides > 6)
        return FillMaskBlock(rowStride, sliceStride, x0, x1, y0, y1, z0, z1, mask);

    if (sides > 0) // -Z
        FillMaskBlock(rowStride, sliceStride, x0, x1, y0, y1, z0, z0, mask);
    if (sides > 1) // -X
        FillMaskBlock(rowStride, sliceStride, x0, x0, y0, y1, z0, z1, mask);
    if (sides > 2) // -Y
        FillMaskBlock(rowStride, sliceStride, x0, x1, y0, y0, z0, z1, mask);
    if (sides > 3) // +X
        FillMaskBlock(rowStride, sliceStride, x1, x1, y0, y1, z0, z1, mask);
    if (sides > 4) // +Y
        FillMaskBlock(rowStride, sliceStride, x0, x1, y1, y1, z0, z1, mask);
    if (sides > 5) // +Z
        FillMaskBlock(rowStride, sliceStride, x0, x1, y0, y1, z1, z1, mask);
}

void DOL::InitDeltasFromBitMask(int w, int h, int d, const uint32_t mask[], cCellDelta3 deltas[], bool invert)
{
    cCellDelta3 maxDelta = { kMaxDelta3, kMaxDelta3, kMaxDelta3 };
    cCellDelta3 minDelta = { 0, 0, 0 };
    
    int s = w * h * d;
    for (int i = 0; i < s; i++)
        deltas[i] = maxDelta;
    
    const int      maskStride = MaskSize(w);
    const uint32_t allClear   = invert ? ~uint32_t(0) : 0;
    
    for (int z = 0; z < d; z++)
        for (int y = 0; y < h; y++)
            for (int j = 0; j < maskStride * 32; j += 32)
            {
                uint32_t m = (*mask++);
                
                if (m == allClear)
                    continue;
                
                cCellDelta3* block = deltas + z * h * w + y * w + j;
                const int n = (w - j) > 32 ? 32 : (w - j);

                for (int i = 0; i < n; i++)
                {
                    if ((m & 1) ^ invert)
                        block[i] = minDelta;
                    
                    m >>= 1;
                }
            }
}

namespace
{
    constexpr cCellDelta3 kB3W0 = { +1, 0, 0 };  // level set is half a cell right
    constexpr cCellDelta3 kB3W1 = { -1, 0, 0 };  // level set is half a cell left
    constexpr cCellDelta3 kB3H0 = { 0, +1, 0 };  // level set is half a cell up
    constexpr cCellDelta3 kB3H1 = { 0, -1, 0 };  // level set is half a cell down
    constexpr cCellDelta3 kB3D0 = { 0, 0, +1 };  // level set is half a cell in
    constexpr cCellDelta3 kB3D1 = { 0, 0, -1 };  // level set is half a cell out

#ifdef USE_BORDER_REFERENCE
    void SetBorderCellsBC(int w, int h, int d, const uint32_t* const mask, cCellDelta3 deltas[])
    {
        int maskStride = MaskSize(w);
        int maskSliceStride = MaskSize(w, h);
        int stride = w;
        int sliceStride = w * h;

        for (int z = 0; z < d - 1; z++)
        {
            const uint32_t* maskSlice0  = mask   + (z + 0) * maskSliceStride;
            const uint32_t* maskSlice1  = mask   + (z + 1) * maskSliceStride;

            cCellDelta3* deltasSlice[2] = { deltas + (z + 0) * sliceStride, deltas + (z + 1) * sliceStride };
        
            const uint32_t* maskRow00   = maskSlice0;
            const uint32_t* maskRow01   = maskSlice0 + maskStride;
            const uint32_t* maskRow10   = maskSlice1;
            const uint32_t* maskRow11   = maskSlice1 + maskStride;

            cCellDelta3* rows[2][2] =
            {
                { deltasSlice[0], deltasSlice[0] + stride },
                { deltasSlice[1], deltasSlice[1] + stride },
            };

            for (int y = 0; y < h - 1; y++)
            {
                cCellDelta3* block[2][2] =
                {
                    { rows[0][0] - 1, rows[0][1] - 1 },
                    { rows[1][0] - 1, rows[1][1] - 1 },
                };

                uint32_t lastBit00 = maskRow00[0] & 1;
                uint32_t lastBit01 = maskRow01[0] & 1;
                uint32_t lastBit10 = maskRow10[0] & 1;
                uint32_t lastBit11 = maskRow11[0] & 1;

                for (int j = 0; j < maskStride; j++)
                {
                    // find all cells in 2x2 box (well, 32 shifted sets of them)
                    uint32_t m000 = (maskRow00[j] << 1) | lastBit00;
                    uint32_t m001 =  maskRow00[j];
                    uint32_t m010 = (maskRow01[j] << 1) | lastBit01;
                    uint32_t m011 =  maskRow01[j];
                    uint32_t m100 = (maskRow10[j] << 1) | lastBit10;
                    uint32_t m101 =  maskRow10[j];
                    uint32_t m110 = (maskRow11[j] << 1) | lastBit11;
                    uint32_t m111 =  maskRow11[j];

                    lastBit00 = m001 >> 31;
                    lastBit01 = m011 >> 31;
                    lastBit10 = m101 >> 31;
                    lastBit11 = m111 >> 31;

                    uint32_t ex =  (m000 ^ m001);
                    uint32_t ey =  (m000 ^ m010);
                    uint32_t ez =  (m000 ^ m100);

                    uint32_t ezx = (m001 ^ m101);
                    uint32_t eyz = (m100 ^ m110);
                    uint32_t exy = (m010 ^ m011);

                    uint32_t fx = (ey ^ eyz); // set if we have a diagonal within x=0 face
                    uint32_t fy = (ez ^ ezx);
                    uint32_t fz = (ex ^ exy);
                    
                    uint32_t fzz = eyz ^ m101 ^ m111; // (m100 ^ m101 ^ m110 ^ m111);
                    uint32_t bc = (fz ^ fzz); 
                    
                    if (bc)
                    for (int i = 0; i < 32; i++)
                    {
                        if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                            break;

                        if (bc & (1 << i))    // necessary but not sufficient
                        {
                            int dx = (fx >> i) & 1;
                            int dy = (fy >> i) & 1;
                            int dz = (fz >> i) & 1;

                            cCellDelta3& corner = block[dz][dy][dx + i];
                            corner.x = 1 - 2 * dx;
                            corner.y = 1 - 2 * dy;
                            corner.z = 1 - 2 * dz;
                        }

                    }

                    block[0][0] += 32;
                    block[0][1] += 32;
                    block[1][0] += 32;
                    block[1][1] += 32;
                }

                maskRow00 += maskStride;
                maskRow01 += maskStride;
                maskRow10 += maskStride;
                maskRow11 += maskStride;

                rows[0][0] += stride;
                rows[0][1] += stride;
                rows[1][0] += stride;
                rows[1][1] += stride;
            }
        }
    }

    void SetBorderCellsFCX(int w, int h, int d, const uint32_t* const mask, cCellDelta3 deltas[])
    {
        int maskStride = MaskSize(w);
        int maskSliceStride = MaskSize(w, h);
        int stride = w;
        int sliceStride = w * h;

        for (int z = 0; z < d - 1; z++)
        {
            const uint32_t* maskSlice0  = mask   + (z + 0) * maskSliceStride;
            const uint32_t* maskSlice1  = mask   + (z + 1) * maskSliceStride;

            cCellDelta3* deltasSlice[2] = { deltas + (z + 0) * sliceStride, deltas + (z + 1) * sliceStride };
        
            const uint32_t* maskRow00   = maskSlice0;
            const uint32_t* maskRow01   = maskSlice0 + maskStride;
            const uint32_t* maskRow10   = maskSlice1;
            const uint32_t* maskRow11   = maskSlice1 + maskStride;

            cCellDelta3* rows[2][2] =
            {
                { deltasSlice[0], deltasSlice[0] + stride },
                { deltasSlice[1], deltasSlice[1] + stride },
            };

            for (int y = 0; y < h - 1; y++)
            {
                cCellDelta3* block[2][2] =
                {
                    { rows[0][0], rows[0][1] },
                    { rows[1][0], rows[1][1] },
                };

                for (int j = 0; j < maskStride; j++)
                {
                    // find all cells in 2x2 box (well, 32 shifted sets of them)
                    uint32_t m000 = maskRow00[j];
                    uint32_t m010 = maskRow01[j];
                    uint32_t m100 = maskRow10[j];
                    uint32_t m110 = maskRow11[j];

                    uint32_t ey =  (m000 ^ m010);
                    uint32_t ez =  (m000 ^ m100);

                    uint32_t eyz = (m100 ^ m110);

                    uint32_t fx = (ey ^ eyz); // set if we have a diagonal within x=0 face

                    if (fx)
                    for (int i = 0; i < 32; i++)
                    {
                        if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                            break;

                        if (fx & (1 << i))
                        {
                            int dy = (ez >> i) & 1;
                            int dz = (ey >> i) & 1;

                            cCellDelta3& corner = block[dz][dy][i];
                            corner.x = 0;
                            corner.y = 1 - 2 * dy;
                            corner.z = 1 - 2 * dz;
                        }
                    }

                    block[0][0] += 32;
                    block[0][1] += 32;
                    block[1][0] += 32;
                    block[1][1] += 32;
                }

                maskRow00 += maskStride;
                maskRow01 += maskStride;
                maskRow10 += maskStride;
                maskRow11 += maskStride;

                rows[0][0] += stride;
                rows[0][1] += stride;
                rows[1][0] += stride;
                rows[1][1] += stride;
            }
        }
    }

    void SetBorderCellsFCY(int w, int h, int d, const uint32_t* const mask, cCellDelta3 deltas[])
    {
        int maskStride = MaskSize(w);
        int maskSliceStride = MaskSize(w, h);
        int stride = w;
        int sliceStride = w * h;

        for (int z = 0; z < d - 1; z++)
        {
            const uint32_t* maskSlice0  = mask   + (z + 0) * maskSliceStride;
            const uint32_t* maskSlice1  = mask   + (z + 1) * maskSliceStride;

            cCellDelta3* deltasSlice[2] = { deltas + (z + 0) * sliceStride, deltas + (z + 1) * sliceStride };
        
            const uint32_t* maskRow00   = maskSlice0;
            const uint32_t* maskRow10   = maskSlice1;

            cCellDelta3* rows[2] =
            {
                deltasSlice[0],
                deltasSlice[1],
            };

            for (int y = 0; y < h; y++)
            {
                cCellDelta3* block[2] =
                {
                    rows[0] - 1,
                    rows[1] - 1,
                };

                uint32_t lastBit00 = maskRow00[0] & 1;
                uint32_t lastBit10 = maskRow10[0] & 1;

                for (int j = 0; j < maskStride; j++)
                {
                    // find all cells in 2x2 box (well, 32 shifted sets of them)
                    uint32_t m000 = (maskRow00[j] << 1) | lastBit00;
                    uint32_t m001 =  maskRow00[j];
                    uint32_t m100 = (maskRow10[j] << 1) | lastBit10;
                    uint32_t m101 =  maskRow10[j];

                    lastBit00 = m001 >> 31;
                    lastBit10 = m101 >> 31;

                    uint32_t ex =  (m000 ^ m001);
                    uint32_t ez =  (m000 ^ m100);
                    uint32_t ezx = (m001 ^ m101);

                    uint32_t fy = (ez ^ ezx);

                    if (fy)
                    for (int i = 0; i < 32; i++)
                    {
                        if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                            break;

                        if (fy & (1 << i))
                        {
                            int dx = (ez >> i) & 1;
                            int dz = (ex >> i) & 1;

                            cCellDelta3& corner = block[dz][dx + i];
                            corner.x = 1 - 2 * dx;
                            corner.y = 0;
                            corner.z = 1 - 2 * dz;
                        }
                    }

                    block[0] += 32;
                    block[1] += 32;
                }

                maskRow00 += maskStride;
                maskRow10 += maskStride;

                rows[0] += stride;
                rows[1] += stride;
            }
        }
    }

    void SetBorderCellsFCZ(int w, int h, int d, const uint32_t* const mask, cCellDelta3 deltas[])
    {
        int maskStride = MaskSize(w);
        int maskSliceStride = MaskSize(w, h);
        int stride = w;
        int sliceStride = w * h;

        for (int z = 0; z < d; z++)
        {
            const uint32_t* maskSlice = mask + z * maskSliceStride;
            cCellDelta3* deltasSlice = deltas + z * sliceStride;
        
            const uint32_t* maskRow00 = maskSlice;
            const uint32_t* maskRow01 = maskSlice + maskStride;

            cCellDelta3* rows[2] = { deltasSlice, deltasSlice + stride };

            for (int y = 0; y < h - 1; y++)
            {
                cCellDelta3* block[2] = { rows[0] - 1, rows[1] - 1 };

                uint32_t lastBit00 = maskRow00[0] & 1;
                uint32_t lastBit01 = maskRow01[0] & 1;

                for (int j = 0; j < maskStride; j++)
                {
                    // find all cells in 2x2 box (well, 32 shifted sets of them)
                    uint32_t m000 = (maskRow00[j] << 1) | lastBit00;
                    uint32_t m001 =  maskRow00[j];
                    uint32_t m010 = (maskRow01[j] << 1) | lastBit01;
                    uint32_t m011 =  maskRow01[j];

                    lastBit00 = m001 >> 31;
                    lastBit01 = m011 >> 31;

                    uint32_t ex =  (m000 ^ m001);
                    uint32_t ey =  (m000 ^ m010);
                    uint32_t exy = (m010 ^ m011);

                    uint32_t fz = (ex ^ exy);
                    
                    if (fz)
                    for (int i = 0; i < 32; i++)
                    {
                        if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                            break;

                        if (fz & (1 << i))
                        {
                            int dx = (ey >> i) & 1;
                            int dy = (ex >> i) & 1;

                            cCellDelta3& corner = block[dy][dx + i];
                            corner.x = 1 - 2 * dx;
                            corner.y = 1 - 2 * dy;
                            corner.z = 0;
                        }
                    }

                    block[0] += 32;
                    block[1] += 32;
                }

                maskRow00 += maskStride;
                maskRow01 += maskStride;

                rows[0] += stride;
                rows[1] += stride;
            }
        }
    }

    void SetBorderCellsD(int w, int h, int d, const uint32_t* const mask, cCellDelta3 deltas[])
    {
        // Do depth edges...
        int maskStride = MaskSize(w);
        int maskSliceStride = MaskSize(w, h);

        for (int z = 0; z < d - 1; z++)
            for (int y = 0; y < h; y++)
            {
                const uint32_t* maskRow0   = mask   + (z + 0) * maskSliceStride + y * maskStride;
                const uint32_t* maskRow1   = mask   + (z + 1) * maskSliceStride + y * maskStride;
                cCellDelta3*    deltasRow0 = deltas + (z + 0) * w * h + y * w;
                cCellDelta3*    deltasRow1 = deltas + (z + 1) * w * h + y * w;

                for (int j = 0; j < maskStride; j++)
                {
                    uint32_t m0 = maskRow0[j];
                    uint32_t m1 = maskRow1[j];
                    uint32_t mx = m0 ^ m1;  // 1 on change

                    if (mx)
                    {
                        cCellDelta3* block0 = deltasRow0 + j * 32;
                        cCellDelta3* block1 = deltasRow1 + j * 32;

                        for (int i = 0; i < 32; i++)
                            if (mx & (1 << i))
                            {
                                if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                                    break;

                                assert(j * 32 + i < w);

                                block0[i] = kB3D0;
                                block1[i] = kB3D1;
                            }
                    }
                }
            }
    }

    void SetBorderCellsH(int w, int h, int d, const uint32_t* const mask, cCellDelta3 deltas[])
    {
        int maskStride = MaskSize(w);
        int maskSliceStride = MaskSize(w, h);

        for (int z = 0; z < d; z++)
        {
            const uint32_t* maskSlice  = mask   + z * maskSliceStride;
            cCellDelta3*    deltaSlice = deltas + z * (w * h);

            // Do horizontal edges...

            for (int y = 0; y < h - 1; y++)
            {
                const uint32_t* maskRow0   = maskSlice  + (y + 0) * maskStride;
                const uint32_t* maskRow1   = maskSlice  + (y + 1) * maskStride;
                cCellDelta3*    deltasRow0 = deltaSlice + (y + 0) * w;
                cCellDelta3*    deltasRow1 = deltaSlice + (y + 1) * w;

                for (int j = 0; j < maskStride; j++)
                {
                    uint32_t m0 = maskRow0[j];
                    uint32_t m1 = maskRow1[j];
                    uint32_t mx = m0 ^ m1;  // 1 on change

                    if (mx)
                    {
                        cCellDelta3* block0 = deltasRow0 + j * 32;
                        cCellDelta3* block1 = deltasRow1 + j * 32;

                        for (int i = 0; i < 32; i++)
                            if (mx & (1 << i))
                            {
                                if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                                    break;

                                assert(j * 32 + i < w);

                                block0[i] = kB3H0;
                                block1[i] = kB3H1;
                            }
                    }
                }
            }
        }
    }
    
    void SetBorderCellsW(int w, int h, int d, const uint32_t* const mask, cCellDelta3 deltas[])
    {
        int maskStride = MaskSize(w);
        int maskSliceStride = MaskSize(w, h);
        
        for (int z = 0; z < d; z++)
        {
            const uint32_t* maskSlice  = mask   + maskSliceStride * z;
            cCellDelta3*    deltaSlice = deltas + w * h * z;

            // Do vertical edges...

            for (int y = 0; y < h; y++)
            {
                const uint32_t* maskRow   = maskSlice  + y * maskStride;
                cCellDelta3*    deltasRow = deltaSlice + y * w;
                
                uint32_t lastBit = maskRow[0] & 1;

                for (int j = 0; j < maskStride; j++)
                {
                    uint32_t m0 = maskRow[j];
                    uint32_t m1 = (m0 << 1) | lastBit;
                    lastBit = m0 >> 31;

                    uint32_t mx = m0 ^ m1;

                    if (mx)
                    {
                        cCellDelta3* block = deltasRow + j * 32;

                        for (int i = 0; i < 32; i++)
                            if (mx & (1 << i))
                            {
                                if (j * 32 + i >= w)    // delayed check as we get here relatively rarely
                                    break;

                                assert(j * 32 + i > 0);
                                assert(j * 32 + i < w);

                                block[i - 1] = kB3W0;
                                block[i + 0] = kB3W1;
                            }
                    }
                }
            }
        }
    }

#else

    void SetBorderCellsCombo(int w, int h, int d, const uint32_t* const mask, cCellDelta3 deltas[])
    {
        // Do depth edges...
        int maskStride = MaskSize(w);
        int maskSliceStride = MaskSize(w, h);
        int stride = w;
        int sliceStride = w * h;

        for (int z = 0; z < d - 1; z++)
        {
            const uint32_t* maskSlice0  = mask   + (z + 0) * maskSliceStride;
            const uint32_t* maskSlice1  = mask   + (z + 1) * maskSliceStride;

            cCellDelta3* deltasSlice[2] = { deltas + (z + 0) * sliceStride, deltas + (z + 1) * sliceStride };
        
            const uint32_t* maskRow00   = maskSlice0;
            const uint32_t* maskRow01   = maskSlice0 + maskStride;
            const uint32_t* maskRow10   = maskSlice1;
            const uint32_t* maskRow11   = maskSlice1 + maskStride;

            cCellDelta3* rows[2][2] =
            {
                { deltasSlice[0], deltasSlice[0] + stride },
                { deltasSlice[1], deltasSlice[1] + stride },
            };

            for (int y = 0; y < h - 1; y++)
            {
                cCellDelta3* block[2][2] =
                {
                    { rows[0][0] - 1, rows[0][1] - 1 },
                    { rows[1][0] - 1, rows[1][1] - 1 },
                };

                uint32_t lastBit00 = maskRow00[0] & 1;
                uint32_t lastBit01 = maskRow01[0] & 1;
                uint32_t lastBit10 = maskRow10[0] & 1;
                uint32_t lastBit11 = maskRow11[0] & 1;

                for (int j = 0; j < maskStride; j++)
                {
                    // find all cells in 2x2 box (well, 32 shifted sets of them)
                    uint32_t m000 = (maskRow00[j] << 1) | lastBit00;
                    uint32_t m001 =  maskRow00[j];
                    uint32_t m010 = (maskRow01[j] << 1) | lastBit01;
                    uint32_t m011 =  maskRow01[j];
                    uint32_t m100 = (maskRow10[j] << 1) | lastBit10;
                    uint32_t m101 =  maskRow10[j];
                    uint32_t m110 = (maskRow11[j] << 1) | lastBit11;
                    uint32_t m111 =  maskRow11[j];

                    lastBit00 = m001 >> 31;
                    lastBit01 = m011 >> 31;
                    lastBit10 = m101 >> 31;
                    lastBit11 = m111 >> 31;

                    uint32_t ex =  (m000 ^ m001);
                    uint32_t ey =  (m000 ^ m010);
                    uint32_t ez =  (m000 ^ m100);

                    uint32_t ezx = (m001 ^ m101);
                    uint32_t eyz = (m100 ^ m110);
                    uint32_t exy = (m010 ^ m011);

                    uint32_t fx = (ey ^ eyz); // set if we have a diagonal within x=0 face
                    uint32_t fy = (ez ^ ezx);
                    uint32_t fz = (ex ^ exy);

                    uint32_t eyx  = (m001 ^ m011);
                    uint32_t eyxz = (m101 ^ m111);
                    
                    uint32_t fxx = eyx ^ eyxz;

                    uint32_t fzz = eyz ^ eyxz;
                    uint32_t bc = (fz ^ fzz);       // Necessary but not sufficient, we're just checking z=0 has three cells the same vs z=1 doesn't (or vice versa). Incorrectly identified box corners will always be overwritten however.  

                    // if (ex | ey | ez | fxx | fy | fz | bc) {   early out makes no appreciable difference
                    {
                        // The order here is important as iterating over the masks destroys them
                        ITER_SET_BITS_BEGIN(bc)
                            assert(j * 32 + i > 0);

                            int dx = (fx >> i) & 1;
                            int dy = (fy >> i) & 1;
                            int dz = (fz >> i) & 1;

                            cCellDelta3& corner = block[dz][dy][dx + i];
                            if (corner.x == kMaxDelta3)
                            {
                                corner.x = 1 - 2 * dx;
                                corner.y = 1 - 2 * dy;
                                corner.z = 1 - 2 * dz;
                            }
                        }

                        ITER_SET_BITS_BEGIN(fxx)
                            int dy = (ezx >> i) & 1;
                            int dz = (eyx >> i) & 1;

                            cCellDelta3& corner = block[dz][dy][1 + i];
                            if (corner.x == kMaxDelta3 || (abs(corner.x) + abs(corner.y) + abs(corner.z) > 1))
                            {
                                corner.x = 0;
                                corner.y = 1 - 2 * dy;
                                corner.z = 1 - 2 * dz;
                            }
                        ITER_SET_BITS_END

                        ITER_SET_BITS_BEGIN(fy)
                            assert(j * 32 + i > 0);

                            int dx = (ez >> i) & 1;
                            int dz = (ex >> i) & 1;

                            cCellDelta3& corner = block[dz][0][dx + i];
                            if (corner.x == kMaxDelta3 || (abs(corner.x) + abs(corner.y) + abs(corner.z) > 1))
                            {
                                corner.x = 1 - 2 * dx;
                                corner.y = 0;
                                corner.z = 1 - 2 * dz;
                            }
                        ITER_SET_BITS_END

                        ITER_SET_BITS_BEGIN(fz)
                            assert(j * 32 + i > 0);

                            int dx = (ey >> i) & 1;
                            int dy = (ex >> i) & 1;

                            cCellDelta3& corner = block[0][dy][dx + i];
                            if (corner.x == kMaxDelta3 || (abs(corner.x) + abs(corner.y) + abs(corner.z) > 1))
                            {
                                corner.x = 1 - 2 * dx;
                                corner.y = 1 - 2 * dy;
                                corner.z = 0;
                            }
                        ITER_SET_BITS_END

                        ITER_SET_BITS_BEGIN(ex)
                            assert(j * 32 + i > 0);

                            block[0][0][i + 0] = kB3W0;
                            block[0][0][i + 1] = kB3W1;
                        ITER_SET_BITS_END

                        ITER_SET_BITS_BEGIN(eyx)
                            block[0][0][i + 1] = kB3H0;
                            block[0][1][i + 1] = kB3H1;
                        ITER_SET_BITS_END

                        ITER_SET_BITS_BEGIN (ezx)
                            block[0][0][i + 1] = kB3D0;
                            block[1][0][i + 1] = kB3D1;
                        ITER_SET_BITS_END
                    }

                    block[0][0] += 32;
                    block[0][1] += 32;
                    block[1][0] += 32;
                    block[1][1] += 32;
                }

                maskRow00 += maskStride;
                maskRow01 += maskStride;
                maskRow10 += maskStride;
                maskRow11 += maskStride;

                rows[0][0] += stride;
                rows[0][1] += stride;
                rows[1][0] += stride;
                rows[1][1] += stride;
            }

            // do last fy without everything else
            {
                cCellDelta3* block[2] =
                {
                    rows[0][0] - 1,
                    rows[1][0] - 1,
                };

                uint32_t lastBit00 = maskRow00[0] & 1;
                uint32_t lastBit10 = maskRow10[0] & 1;

                for (int j = 0; j < maskStride; j++)
                {
                    uint32_t m000 = (maskRow00[j] << 1) | lastBit00;
                    uint32_t m001 =  maskRow00[j];
                    uint32_t m100 = (maskRow10[j] << 1) | lastBit10;
                    uint32_t m101 =  maskRow10[j];

                    lastBit00 = m001 >> 31;
                    lastBit10 = m101 >> 31;

                    uint32_t ex =  (m000 ^ m001);
                    uint32_t ez =  (m000 ^ m100);
                    uint32_t ezx = (m001 ^ m101);

                    uint32_t fy = (ez ^ ezx);
                    
                    ITER_SET_BITS_BEGIN(fy)
                        assert(j * 32 + i > 0);

                        int dx = (ez >> i) & 1;
                        int dz = (ex >> i) & 1;

                        cCellDelta3& corner = block[dz][dx + i];
                        if (corner.x == kMaxDelta3 || (abs(corner.x) + abs(corner.y) + abs(corner.z) > 1))
                        {
                            corner.x = 1 - 2 * dx;
                            corner.y = 0;
                            corner.z = 1 - 2 * dz;
                        }
                    ITER_SET_BITS_END

                    ITER_SET_BITS_BEGIN(ex)
                        assert(j * 32 + i > 0);

                        block[0][i + 0] = kB3W0;
                        block[0][i + 1] = kB3W1;
                    ITER_SET_BITS_END

                    ITER_SET_BITS_BEGIN(ezx)
                        block[0][i + 1] = kB3D0;
                        block[1][i + 1] = kB3D1;
                    ITER_SET_BITS_END

                    block[0] += 32;
                    block[1] += 32;
                }
            }
        }
        
        // Now for the last slice, do the fz faces, and ex/ey edges
        {
            int z = d - 1;
            const uint32_t* maskSlice0  = mask   + z * maskSliceStride;
            cCellDelta3*    deltasSlice = deltas + z * sliceStride;
        
            const uint32_t* maskRow00   = maskSlice0;
            const uint32_t* maskRow01   = maskSlice0 + maskStride;

            cCellDelta3* rows[2] = { deltasSlice, deltasSlice + stride };

            for (int y = 0; y < h - 1; y++)
            {
                cCellDelta3* block[2] = { rows[0] - 1, rows[1] - 1 };

                uint32_t lastBit00 = maskRow00[0] & 1;
                uint32_t lastBit01 = maskRow01[0] & 1;

                for (int j = 0; j < maskStride; j++)
                {
                    uint32_t m000 = (maskRow00[j] << 1) | lastBit00;
                    uint32_t m001 =  maskRow00[j];
                    uint32_t m010 = (maskRow01[j] << 1) | lastBit01;
                    uint32_t m011 =  maskRow01[j];

                    lastBit00 = m001 >> 31;
                    lastBit01 = m011 >> 31;

                    uint32_t ex =  (m000 ^ m001);
                    uint32_t ey =  (m000 ^ m010);
                    uint32_t eyx = (m001 ^ m011);

                    uint32_t fz = ey ^ eyx;

                    ITER_SET_BITS_BEGIN(fz)
                        assert(j * 32 + i > 0);

                        int dx = (ey >> i) & 1;
                        int dy = (ex >> i) & 1;

                        assert(i != 0 || j != 0 || dx != 0);

                        cCellDelta3& corner = block[dy][dx + i];
                        if (corner.x == kMaxDelta3 || (abs(corner.x) + abs(corner.y) + abs(corner.z) > 1))
                        {
                            corner.x = 1 - 2 * dx;
                            corner.y = 1 - 2 * dy;
                            corner.z = 0;
                        }
                    ITER_SET_BITS_END

                    ITER_SET_BITS_BEGIN(ex)
                        assert(j * 32 + i > 0);

                        block[0][i + 0] = kB3W0;
                        block[0][i + 1] = kB3W1;
                    ITER_SET_BITS_END

                    ITER_SET_BITS_BEGIN(eyx)
                        block[0][i + 1] = kB3H0;
                        block[1][i + 1] = kB3H1;
                    ITER_SET_BITS_END

                    block[0] += 32;
                    block[1] += 32;
                }

                maskRow00 += maskStride;
                maskRow01 += maskStride;

                rows[0] += stride;
                rows[1] += stride;
            }
            
            // do last ex without everything else
            cCellDelta3* block = rows[0] - 1;
            uint32_t lastBit00 = maskRow00[0] & 1;

            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m000 = (maskRow00[j] << 1) | lastBit00;
                uint32_t m001 =  maskRow00[j];

                lastBit00 = m001 >> 31;

                uint32_t ex =  (m000 ^ m001);

                ITER_SET_BITS_BEGIN(ex)
                    assert(j * 32 + i > 0);

                    block[i + 0] = kB3W0;
                    block[i + 1] = kB3W1;
                ITER_SET_BITS_END

                block += 32;
            }
        }
    }
#endif
}

void DOL::InitBorderDeltasFromBitMask(int w, int h, int d, const uint32_t mask[], cCellDelta3 deltas[])
{
    const cCellDelta3 maxDelta = { kMaxDelta3, kMaxDelta3, kMaxDelta3 };

    int s = w * h * d;
    for (int i = 0; i < s; i++)
        deltas[i] = maxDelta;

#ifdef USE_BORDER_REFERENCE
    SetBorderCellsBC (w, h, d, mask, deltas);
    SetBorderCellsFCZ(w, h, d, mask, deltas);
    SetBorderCellsFCY(w, h, d, mask, deltas);
    SetBorderCellsFCX(w, h, d, mask, deltas);
    SetBorderCellsD  (w, h, d, mask, deltas);
    SetBorderCellsH  (w, h, d, mask, deltas);
    SetBorderCellsW  (w, h, d, mask, deltas);
#else
    SetBorderCellsCombo(w, h, d, mask, deltas);
#endif
}

namespace
{
    inline bool Closer(cCellDelta3 d0, cCellDelta3 d1)
    {
        int dw0 = d0.x * d0.x + d0.y * d0.y + d0.z * d0.z;
        int dw1 = d1.x * d1.x + d1.y * d1.y + d1.z * d1.z;
        
        return dw0 > dw1;
    }
}

void DOL::FastSweep(int w, int h, int d, cCellDelta3 deltas[], int cw)
{
    // Fast sweep approach -- sweep each diagonal in turn.
    // C B
    // A x
    const int sliceStride = w * h;
    const int lastRow = w * (h - 1);

    for (int sweep = 0; sweep < 8; sweep++)
    {
        const int sx = (sweep >> 0) & 1;
        const int sy = (sweep >> 1) & 1;
        const int sz = (sweep >> 2) & 1;

        const int dx = 1 - 2 * sx;
        const int dy = 1 - 2 * sy;
        const int dz = 1 - 2 * sz;

        const int cx = -cw * dx;
        const int cy = -cw * dy;
        const int cz = -cw * dz;

        const int ib =     sx * (w - 1);
        const int ie = w - sx * (w + 1);

        cCellDelta3* slice0 = deltas + sz * (d - 1) * sliceStride;
        cCellDelta3 cell;

        for (int k = 1; k < d; k++)
        {
            cCellDelta3* slice1 = slice0 + sliceStride * dz;
            
            cCellDelta3* row00 = slice0 + sy * lastRow;
            cCellDelta3* row01 = slice1 + sy * lastRow;
            
            for (int j = 1; j < h; j++)
            {
                cCellDelta3* row10 = row00 + w * dy;
                cCellDelta3* row11 = row01 + w * dy;
            
                for (int i = ib; i != ie; i += dx)
                {
                    cell.x = row10[i].x;
                    cell.y = row10[i].y;
                    cell.z = row10[i].z + cz;
                    
                    if (Closer(row11[i], cell))
                        row11[i] = cell;

                    cell.x = row01[i].x;
                    cell.y = row01[i].y + cy;
                    cell.z = row01[i].z;
                    
                    if (Closer(row11[i], cell))
                        row11[i] = cell;
                }

                for (int i = ib + dx; i != ie; i += dx)
                {
                    cell.x = row11[i - dx].x + cx;
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

namespace
{
    void DanielssonSlice(int w, int h, cCellDelta3 deltas[], int cw)
    {
        cCellDelta3* row0 = deltas;
        cCellDelta3 cell;

        for (int j = 1; j < h; j++)
        {
            cCellDelta3* row1 = row0 + w;

            for (int i = 0; i < w; i++)
            {
                cell.x = row0[i].x;
                cell.y = row0[i].y - cw;
                cell.z = row0[i].z;

                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }

            for (int i = 1; i < w; i++)
            {
                cell.x = row1[i - 1].x - cw;
                cell.y = row1[i - 1].y;
                cell.z = row1[i - 1].z;

                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }

            for (int i = w - 2; i >= 0; i--)
            {

                cell.x = row1[i + 1].x + cw;
                cell.y = row1[i + 1].y;
                cell.z = row1[i + 1].z;

                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }

            row0 = row1;
        }

        row0 = deltas + h * w - w;

        for (int j = h - 1; j > 0; j--)
        {
            cCellDelta3* row1 = row0 - w;

            for (int i = 0; i < w; i++)
            {
                cell.x = row0[i].x;
                cell.y = row0[i].y + cw;
                cell.z = row0[i].z;

                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }

            for (int i = 1; i < w; i++)
            {
                cell.x = row1[i - 1].x - cw;
                cell.y = row1[i - 1].y;
                cell.z = row1[i - 1].z;

                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }

            for (int i = w - 2; i >= 0; i--)
            {
                cell.x = row1[i + 1].x + cw;
                cell.y = row1[i + 1].y;
                cell.z = row1[i + 1].z;

                if (Closer(row1[i], cell))
                    row1[i] = cell;
            }

            row0 = row1;
        }
    }
}

void DOL::Danielsson(int w, int h, int d, cCellDelta3 deltas[], int cw)
{
    int sliceStride = w * h;

    cCellDelta3* slice0 = deltas;
    cCellDelta3 cell;

    for (int k = 1; k < d; k++)
    {
        cCellDelta3* slice1 = slice0 + sliceStride;

        for (int i = 0; i < sliceStride; i++)
        {
            cell.x = slice0[i].x;
            cell.y = slice0[i].y;
            cell.z = slice0[i].z - cw;

            if (Closer(slice1[i], cell))
                slice1[i] = cell;
        }

        DanielssonSlice(w, h, slice1, cw);

        slice0 = slice1;
    }

    slice0 = deltas + sliceStride * (d - 1);

    for (int k = d - 2; k >= 0; k--)
    {
        cCellDelta3* slice1 = slice0 - sliceStride;

        for (int i = 0; i < sliceStride; i++)
        {
            cell.x = slice0[i].x;
            cell.y = slice0[i].y;
            cell.z = slice0[i].z + cw;

            if (Closer(slice1[i], cell))
                slice1[i] = cell;
        }

        DanielssonSlice(w, h, slice1, cw);

        slice0 = slice1;
    }
}

namespace
{
    void JumpFlood(int step, int w, int h, int d, const cCellDelta3 deltasIn[], cCellDelta3 deltasOut[], int cw)
    {
        const cCellDelta3* row    = deltasIn;
        cCellDelta3*       rowOut = deltasOut;

        int wh = w * h;

        for (int k = 0; k < d; k++)
        {
            int kc0 = (k - step) >= 0 ? -step : 0;
            int kc1 = (k + step) <  d ? +step : 0;

            for (int j = 0; j < h; j++)
            {
                int jc0 = (j - step) >= 0 ? -step : 0;
                int jc1 = (j + step) <  h ? +step : 0;

                for (int i = 0; i < w; i++)
                {
                    cCellDelta3 minCell = row[i];

                    int ic0 = (i - step) >= 0 ? -step : 0;
                    int ic1 = (i + step) <  w ? +step : 0;

                    for (int kc = kc0; kc <= kc1; kc += step)
                    for (int jc = jc0; jc <= jc1; jc += step)
                    for (int ic = ic0; ic <= ic1; ic += step)
                    {
                        cCellDelta3 cell = row[i + ic + jc * w + kc * wh];
                        
                        if (cell.x != kMaxDelta3)
                        {
                            cell.x += cw * ic;
                            cell.y += cw * jc;
                            cell.z += cw * kc;

                            if (Closer(minCell, cell))
                                minCell = cell;
                        }
                    }

                    rowOut[i] = minCell;
                }
                
                row    += w;
                rowOut += w;
            }
        }
    }
}

void DOL::JumpFlood(int w, int h, int d, cCellDelta3 deltas[], int cw)
{
    cCellDelta3* temp = new cCellDelta3[w * h * d];
    cCellDelta3* d0 = deltas;
    cCellDelta3* d1 = temp;

    int wh = w >= h ? w : h;
    int step = FloorPow2((wh >= d ? wh : d) - 1);
    
    for ( ; step > 0; step /= 2)
    {
        ::JumpFlood(step, w, h, d, d0, d1, cw);

        cCellDelta3* dt = d0;
        d0 = d1;
        d1 = dt;
    }

    if (d0 != deltas)
        ::JumpFlood(1, w, h, d, d0, deltas, cw);    // we could copy, but might as well do another round

    delete[] temp;
}


// Brute force -- for checking only!
void DOL::BruteForce(int w, int h, int d, cCellDelta3 deltas[], int cw)
{
    for (int k = 0; k < d; k++)
        for (int j = 0; j < h; j++)
            for (int i = 0; i < w; i++)
            {
                cCellDelta3* cell = deltas + k * w * h + j * w + i;
                
                if (abs(cell->x) >= cw || abs(cell->y) >= cw || abs(cell->z) >= cw) // non-border cell
                    continue;

                // For each occupied cell, iterate over all other cells and check for minimal distance
                for (int z = 0; z < d; z++)
                    for (int y = 0; y < h; y++)
                        for (int x = 0; x < w; x++)
                        {
                            cCellDelta3* cell = deltas + z * w * h + y * w + x;
                            
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
            
        const int ib =     sx * (w - 1);
        const int ie = w - sx * (w + 1);

        const int wh  = w * h;      // slice stride
        const int whd = w * h * d;  // volume stride

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
            
        const int ib =     sx * (w - 1);
        const int ie = w - sx * (w + 1);

        const int wh  = w * h;
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
