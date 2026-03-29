//
// DistanceFields.cpp
//
// Distance field generation utilities
//
// Andrew Willmott
//

#include "DistanceFields.hpp"

#if defined(DOS_DEBUG) || defined(DEBUG)
    #include <assert.h>
    #define DOL_ASSERT assert
#else
    #define DOL_ASSERT(e) ((void) 0)
#endif

#if !defined(__cplusplus) || (__cplusplus <= 201100)
    #define constexpr const
#endif

using namespace DOL;

// Declarations

// Convenience macro for iterating over set bits (and also testing different strategies)
#ifdef _MSC_VER
    #include <intrin.h>
    uint32_t TrailingZeroes32(uint32_t mask)
    {
        unsigned long index = 32;
        _BitScanForward(&index, mask);  // index of lowest 1 bit. doesn't set for mask=0, hence default of 32 above
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

#define ITER_SET_BITS_END }

namespace
{
    template<class T> void InitFromBitMask(int w, int h, const uint32_t mask[], T data[], T minD, T maxD)
    {
        // mask is expected to have rows padded to whole numbers of uint32_ts
        int s = w * h;
        for (int i = 0; i < s; i++)
            data[i] = maxD;

        int maskStride = MaskSize(w);

        for (int y = 0; y < h; y++)
        {
            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m = *mask++;

                if (m == 0)
                    continue;

                T* block = data + y * w + 32 * j;

                ITER_SET_BITS_BEGIN(m)
                    block[i] = minD;
                ITER_SET_BITS_END
            }
        }
    }
}

// --- 2D Images ---------------------------------------------------------------

void DOL::InitDistancesFromBitMask(int w, int h, const uint32_t mask[], int16_t distances[], int16_t maxD)
{
    return ::InitFromBitMask<int16_t>(w, h, mask, distances, 0, maxD);
}

void DOL::InitDistancesFromBitMask(int w, int h, const uint32_t mask[], int32_t distances[], int32_t maxD)
{
    return ::InitFromBitMask<int32_t>(w, h, mask, distances, 0, maxD);
}

void DOL::InitDistancesFromBitMask(int w, int h, const uint32_t mask[], float distances[], float maxD)
{
    return ::InitFromBitMask<float>(w, h, mask, distances, 0.0f, maxD);
}

void DOL::InitDistancesFromBitMask(int w, int h, const uint32_t mask[], float distances[], float minD, float maxD)
{
    InitFromBitMask<float>(w, h, mask, distances, minD, maxD);
}

void DOL::InitDeltasFromBitMask(int w, int h, const uint32_t mask[], CellDelta2 deltas[], bool invert)
{
    const CellDelta2 maxDelta = { kMaxDelta2, kMaxDelta2 };
    const CellDelta2 minDelta = { 0, 0 };

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

            CellDelta2* block = deltas + y * w + j;
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
    // the border. (The standard technique of calculating the exterior distances,
    // then inverting the problem to calculate interior distances, and combining
    // the two results, is both 2x as expensive, and leads to inconsistent
    // distances on border cells. The border effectively has a width of -ve 1.

    constexpr CellDelta2 kBH0 = { 0, +1 };  // level set is half a cell up
    constexpr CellDelta2 kBH1 = { 0, -1 };  // level set is half a cell down
    constexpr CellDelta2 kBW0 = { +1, 0 };  // level set is half a cell right
    constexpr CellDelta2 kBW1 = { -1, 0 };  // level set is half a cell left

    void SetBorderCellsCombo(int w, int h, const uint32_t* const mask, CellDelta2 deltas[])
    {
        // Rolls everything into one pass
        const int maskStride = MaskSize(w);

        for (int y = 0; y < h - 1; y++)
        {
            const uint32_t* maskRow0   = mask   + (y + 0) * maskStride;
            const uint32_t* maskRow1   = mask   + (y + 1) * maskStride;
            CellDelta2*     deltasRow0 = deltas + (y + 0) * w;
            CellDelta2*     deltasRow1 = deltas + (y + 1) * w;

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

                CellDelta2* block[2] = { deltasRow0 + j * 32, deltasRow1 + j * 32 };

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
                    DOL_ASSERT(j * 32 + i > 0);

                    int dx = (eyx >> i) & 1;
                    int dy = (ex  >> i) & 1;

                    CellDelta2& corner = block[dy][i + dx - 1];

                    if (corner.x != kMaxDelta2)
                        continue;

                    corner.x = 1 - 2 * dx;
                    corner.y = 1 - 2 * dy;
                ITER_SET_BITS_END

                ITER_SET_BITS_BEGIN(ex)
                    DOL_ASSERT(j * 32 + i > 0);
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
            CellDelta2*     deltasRow = deltas + y * w;

            uint32_t lastBit = maskRow[0] & 1; // repeat border cell, avoids extra logic to loop from 1 rather than 0.

            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m0 = maskRow[j];
                uint32_t m1 = (m0 << 1) | lastBit;
                lastBit = m0 >> 31;

                uint32_t mvx = m0 ^ m1;

                if (mvx)
                {
                    CellDelta2* block = deltasRow + j * 32;

                    ITER_SET_BITS_BEGIN(mvx)
                        DOL_ASSERT(j * 32 + i > 0);
                        block[i - 1] = kBW0;
                        block[i + 0] = kBW1;
                    ITER_SET_BITS_END
                }
            }
        }
    }
}

void DOL::InitBorderDeltasFromBitMask(int w, int h, const uint32_t* const mask, CellDelta2 deltas[])
{
    const CellDelta2 maxDelta = { kMaxDelta2, kMaxDelta2 };

    int s = w * h;
    for (int i = 0; i < s; i++)
        deltas[i] = maxDelta;

    SetBorderCellsCombo(w, h, mask, deltas);
}

namespace
{
    /// Returns true if d1 is closer than d0
    inline bool Closer(CellDelta2 d0, CellDelta2 d1)
    {
        int dw0 = d0.x * d0.x + d0.y * d0.y;
        int dw1 = d1.x * d1.x + d1.y * d1.y;

        return dw0 > dw1;
    }
}

void DOL::FastSweep(int w, int h, CellDelta2 deltas[], int cw)
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
    CellDelta2 cell;

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

        CellDelta2* row0 = deltas + sy * lastRow;

        for (int j = 1; j < h; j++)
        {
            CellDelta2* row1 = row0 + w * dy;

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

void DOL::Danielsson(int w, int h, CellDelta2 deltas[], int cw)
{
    // Store dx^2, dy^2 to nearest point. Distance = sqrt(dx^2 + dy^2).
    // Much more accurate, doesn't over-estimate like Chamfer.
    CellDelta2* row0 = deltas;
    CellDelta2 cell;

    for (int j = 1; j < h; j++)
    {
        CellDelta2* row1 = row0 + w;

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
        CellDelta2* row1 = row0 - w;

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
    void JumpFlood(int step, int w, int h, CellDelta2 deltasIn[], CellDelta2 deltasOut[], int cw)
    {
        #pragma omp parallel for
        for (int j = 0; j < h; j++)
        {
            const CellDelta2* row    = deltasIn  + j * w;
            CellDelta2*       rowOut = deltasOut + j * w;

            for (int i = 0; i < w; i++)
            {
                CellDelta2 minCell = row[i];

                if (i - step >= 0)
                {
                    CellDelta2 cell = row[i - step];
                    cell.x -= cw * step;

                    if (Closer(minCell, cell))
                        minCell = cell;

                    if (j - step >= 0)
                    {
                        CellDelta2 cell = row[i - step - step * w];
                        cell.x -= cw * step;
                        cell.y -= cw * step;

                        if (Closer(minCell, cell))
                            minCell = cell;
                    }

                    if (j + step < h)
                    {
                        CellDelta2 cell = row[i - step + step * w];
                        cell.x -= cw * step;
                        cell.y += cw * step;

                        if (Closer(minCell, cell))
                            minCell = cell;
                    }
                }

                if (i + step < w)
                {
                    CellDelta2 cell = row[i + step];
                    cell.x += cw * step;

                    if (Closer(minCell, cell))
                        minCell = cell;

                    if (j - step >= 0)
                    {
                        CellDelta2 cell = row[i + step - step * w];
                        cell.x += cw * step;
                        cell.y -= cw * step;

                        if (Closer(minCell, cell))
                            minCell = cell;
                    }

                    if (j + step < h)
                    {
                        CellDelta2 cell = row[i + step + step * w];
                        cell.x += cw * step;
                        cell.y += cw * step;

                        if (Closer(minCell, cell))
                            minCell = cell;
                    }
                }

                if (j - step >= 0)
                {
                    CellDelta2 cell = row[i - step * w];
                    cell.y -= cw * step;

                    if (Closer(minCell, cell))
                        minCell = cell;
                }

                if (j + step < h)
                {
                    CellDelta2 cell = row[i + step * w];
                    cell.y += cw * step;

                    if (Closer(minCell, cell))
                        minCell = cell;
                }

                rowOut[i] = minCell;
            }
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

void DOL::JumpFlood(int w, int h, CellDelta2 deltas[], int cw)
{
    CellDelta2* temp = new CellDelta2[w * h];
    CellDelta2* d0 = deltas;
    CellDelta2* d1 = temp;

    int step = ::FloorPow2((w >= h ? w : h) - 1);

    for ( ; step > 0; step /= 2)
    {
        ::JumpFlood(step, w, h, d0, d1, cw);

        CellDelta2* dt = d0;
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
    void Felzenszwalb(int n, int m, const float f2[/* i . m */], float d[], int v[/* n */], float z[/* n + 1 */])
    {
        // See algorithm 1., same variable names are used.
        int k = 0;

        v[0] = 0;
        z[0] = -FLT_MAX;
        z[1] = +FLT_MAX;

        for (int q = 1; q < n; q++)
        {
            float fqt = f2[m * q] + q * q;
            float fvt = f2[m * v[k]] + (v[k] * v[k]);
            float s = (fqt - fvt) / (2 * (q - v[k]));

            while (s <= z[k])
            {
                k--;
                fvt = f2[m * v[k]] + (v[k] * v[k]);
                s = (fqt - fvt) / (2 * (q - v[k]));
            }

            k++;
            v[k] = q;
            z[k] = s;
            z[k + 1] = +FLT_MAX;
        }

        k = 0;
        for (int q = 0; q < n; q++)
        {
            while (z[k + 1] < q)
                k++;

            d[q] = (q - v[k]) * (q - v[k]) + f2[m * v[k]];
        }
    }
}

void DOL::Felzenszwalb(int w, int h, float distances[])
{
    int whMax = w >= h ? w : h;
    float* td = new float[w * h];

    // Transform along columns
    #pragma omp parallel
    {
        int*   v = new int  [whMax];
        float* z = new float[whMax + 1];

        #pragma omp for
        for (int x = 0; x < w; x++)
            ::Felzenszwalb(h, w, distances + x, td + x * h, v, z);

        #pragma omp for
        for (int y = 0; y < h; y++)
            ::Felzenszwalb(w, h, td + y, distances + y * w, v, z);

        delete[] v;
        delete[] z;
    }

    delete[] td;
}

// Brute force -- for checking only!
void DOL::BruteForce(int w, int h, CellDelta2 deltas[], int cw)
{
    for (int i = 0; i < h; i++)
    for (int j = 0; j < w; j++)
    {
        CellDelta2* cell = deltas + i * w + j;

        if (abs(cell->x) >= cw || abs(cell->y) >= cw)   // non-boundary cell?
            continue;

        // For each occupied cell, iterate over all other cells and check for minimal distance

        for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            CellDelta2 d2 = { int16_t((j - x) * cw + cell->x), int16_t((i - y) * cw + cell->y) };
            CellDelta2* cell2 = deltas + y * w + x;

            if (Closer(*cell2, d2))
                *cell2 = d2;
        }
    }
}

// --- 3D Volumes --------------------------------------------------------------

void DOL::InitDeltasFromBitMask(int w, int h, int d, const uint32_t mask[], CellDelta3 deltas[], bool invert)
{
    CellDelta3 maxDelta = { kMaxDelta3, kMaxDelta3, kMaxDelta3 };
    CellDelta3 minDelta = { 0, 0, 0 };

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

                CellDelta3* block = deltas + z * h * w + y * w + j;
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
    constexpr CellDelta3 kB3W0 = { +1, 0, 0 };  // level set is half a cell right
    constexpr CellDelta3 kB3W1 = { -1, 0, 0 };  // level set is half a cell left
    constexpr CellDelta3 kB3H0 = { 0, +1, 0 };  // level set is half a cell up
    constexpr CellDelta3 kB3H1 = { 0, -1, 0 };  // level set is half a cell down
    constexpr CellDelta3 kB3D0 = { 0, 0, +1 };  // level set is half a cell in
    constexpr CellDelta3 kB3D1 = { 0, 0, -1 };  // level set is half a cell out

    void SetBorderCellsCombo(int w, int h, int d, const uint32_t* const mask, CellDelta3 deltas[])
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

            CellDelta3* deltasSlice[2] = { deltas + (z + 0) * sliceStride, deltas + (z + 1) * sliceStride };

            const uint32_t* maskRow00   = maskSlice0;
            const uint32_t* maskRow01   = maskSlice0 + maskStride;
            const uint32_t* maskRow10   = maskSlice1;
            const uint32_t* maskRow11   = maskSlice1 + maskStride;

            CellDelta3* rows[2][2] =
            {
                { deltasSlice[0], deltasSlice[0] + stride },
                { deltasSlice[1], deltasSlice[1] + stride },
            };

            for (int y = 0; y < h - 1; y++)
            {
                CellDelta3* block[2][2] =
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

                    // if (ex | ey | ez | fxx | fy | fz | bc)  //  early out makes no appreciable difference
                    {
                        // The order here is important as iterating over the masks destroys them
                        ITER_SET_BITS_BEGIN(bc)
                            DOL_ASSERT(j * 32 + i > 0);

                            int dx = (fx >> i) & 1;
                            int dy = (fy >> i) & 1;
                            int dz = (fz >> i) & 1;

                            CellDelta3& corner = block[dz][dy][dx + i];
                            if (corner.x == kMaxDelta3)
                            {
                                corner.x = 1 - 2 * dx;
                                corner.y = 1 - 2 * dy;
                                corner.z = 1 - 2 * dz;
                            }
                        ITER_SET_BITS_END

                        ITER_SET_BITS_BEGIN(fxx)
                            int dy = (ezx >> i) & 1;
                            int dz = (eyx >> i) & 1;

                            CellDelta3& corner = block[dz][dy][1 + i];
                            if (corner.x == kMaxDelta3 || (abs(corner.x) + abs(corner.y) + abs(corner.z) > 1))
                            {
                                corner.x = 0;
                                corner.y = 1 - 2 * dy;
                                corner.z = 1 - 2 * dz;
                            }
                        ITER_SET_BITS_END

                        ITER_SET_BITS_BEGIN(fy)
                            DOL_ASSERT(j * 32 + i > 0);

                            int dx = (ez >> i) & 1;
                            int dz = (ex >> i) & 1;

                            CellDelta3& corner = block[dz][0][dx + i];
                            if (corner.x == kMaxDelta3 || (abs(corner.x) + abs(corner.y) + abs(corner.z) > 1))
                            {
                                corner.x = 1 - 2 * dx;
                                corner.y = 0;
                                corner.z = 1 - 2 * dz;
                            }
                        ITER_SET_BITS_END

                        ITER_SET_BITS_BEGIN(fz)
                            DOL_ASSERT(j * 32 + i > 0);

                            int dx = (ey >> i) & 1;
                            int dy = (ex >> i) & 1;

                            CellDelta3& corner = block[0][dy][dx + i];
                            if (corner.x == kMaxDelta3 || (abs(corner.x) + abs(corner.y) + abs(corner.z) > 1))
                            {
                                corner.x = 1 - 2 * dx;
                                corner.y = 1 - 2 * dy;
                                corner.z = 0;
                            }
                        ITER_SET_BITS_END

                        ITER_SET_BITS_BEGIN(ex)
                            DOL_ASSERT(j * 32 + i > 0);

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
                CellDelta3* block[2] =
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
                        DOL_ASSERT(j * 32 + i > 0);

                        int dx = (ez >> i) & 1;
                        int dz = (ex >> i) & 1;

                        CellDelta3& corner = block[dz][dx + i];
                        if (corner.x == kMaxDelta3 || (abs(corner.x) + abs(corner.y) + abs(corner.z) > 1))
                        {
                            corner.x = 1 - 2 * dx;
                            corner.y = 0;
                            corner.z = 1 - 2 * dz;
                        }
                    ITER_SET_BITS_END

                    ITER_SET_BITS_BEGIN(ex)
                        DOL_ASSERT(j * 32 + i > 0);

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
            CellDelta3*    deltasSlice = deltas + z * sliceStride;

            const uint32_t* maskRow00   = maskSlice0;
            const uint32_t* maskRow01   = maskSlice0 + maskStride;

            CellDelta3* rows[2] = { deltasSlice, deltasSlice + stride };

            for (int y = 0; y < h - 1; y++)
            {
                CellDelta3* block[2] = { rows[0] - 1, rows[1] - 1 };

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
                        DOL_ASSERT(j * 32 + i > 0);

                        int dx = (ey >> i) & 1;
                        int dy = (ex >> i) & 1;

                        DOL_ASSERT(i != 0 || j != 0 || dx != 0);

                        CellDelta3& corner = block[dy][dx + i];
                        if (corner.x == kMaxDelta3 || (abs(corner.x) + abs(corner.y) + abs(corner.z) > 1))
                        {
                            corner.x = 1 - 2 * dx;
                            corner.y = 1 - 2 * dy;
                            corner.z = 0;
                        }
                    ITER_SET_BITS_END

                    ITER_SET_BITS_BEGIN(ex)
                        DOL_ASSERT(j * 32 + i > 0);

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
            CellDelta3* block = rows[0] - 1;
            uint32_t lastBit00 = maskRow00[0] & 1;

            for (int j = 0; j < maskStride; j++)
            {
                uint32_t m000 = (maskRow00[j] << 1) | lastBit00;
                uint32_t m001 =  maskRow00[j];

                lastBit00 = m001 >> 31;

                uint32_t ex =  (m000 ^ m001);

                ITER_SET_BITS_BEGIN(ex)
                    DOL_ASSERT(j * 32 + i > 0);

                    block[i + 0] = kB3W0;
                    block[i + 1] = kB3W1;
                ITER_SET_BITS_END

                block += 32;
            }
        }
    }
}

void DOL::InitBorderDeltasFromBitMask(int w, int h, int d, const uint32_t mask[], CellDelta3 deltas[])
{
    const CellDelta3 maxDelta = { kMaxDelta3, kMaxDelta3, kMaxDelta3 };

    int s = w * h * d;
    for (int i = 0; i < s; i++)
        deltas[i] = maxDelta;

    SetBorderCellsCombo(w, h, d, mask, deltas);
}

namespace
{
    inline bool Closer(CellDelta3 d0, CellDelta3 d1)
    {
        int dw0 = d0.x * d0.x + d0.y * d0.y + d0.z * d0.z;
        int dw1 = d1.x * d1.x + d1.y * d1.y + d1.z * d1.z;

        return dw0 > dw1;
    }
}

void DOL::FastSweep(int w, int h, int d, CellDelta3 deltas[], int cw)
{
    // Fast sweep approach -- sweep each diagonal in turn.
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

        CellDelta3* slice0 = deltas + sz * (d - 1) * sliceStride;
        CellDelta3 cell;

        for (int k = 1; k < d; k++)
        {
            CellDelta3* slice1 = slice0 + sliceStride * dz;

            CellDelta3* row00 = slice0 + sy * lastRow;
            CellDelta3* row01 = slice1 + sy * lastRow;

            for (int j = 1; j < h; j++)
            {
                CellDelta3* row10 = row00 + w * dy;
                CellDelta3* row11 = row01 + w * dy;

                for (int i = ib; i != ie; i += dx)
                {
                    // z-neighbor (prev slice, same row)
                    cell.x = row10[i].x;
                    cell.y = row10[i].y;
                    cell.z = row10[i].z + cz;

                    if (Closer(row11[i], cell))
                        row11[i] = cell;

                    // y-neighbor (same slice, prev row)
                    cell.x = row01[i].x;
                    cell.y = row01[i].y + cy;
                    cell.z = row01[i].z;

                    if (Closer(row11[i], cell))
                        row11[i] = cell;

                    // yz-diagonal (prev slice, prev row)
                    cell.x = row00[i].x;
                    cell.y = row00[i].y + cy;
                    cell.z = row00[i].z + cz;

                    if (Closer(row11[i], cell))
                        row11[i] = cell;
                }

                for (int i = ib + dx; i != ie; i += dx)
                {
                    // x-neighbor (same slice, same row, prev col)
                    cell.x = row11[i - dx].x + cx;
                    cell.y = row11[i - dx].y;
                    cell.z = row11[i - dx].z;

                    if (Closer(row11[i], cell))
                        row11[i] = cell;

                    // xz-diagonal (prev slice, same row, prev col)
                    cell.x = row10[i - dx].x + cx;
                    cell.y = row10[i - dx].y;
                    cell.z = row10[i - dx].z + cz;

                    if (Closer(row11[i], cell))
                        row11[i] = cell;

                    // xy-diagonal (same slice, prev row, prev col)
                    cell.x = row01[i - dx].x + cx;
                    cell.y = row01[i - dx].y + cy;
                    cell.z = row01[i - dx].z;

                    if (Closer(row11[i], cell))
                        row11[i] = cell;

                    // xyz-diagonal (prev slice, prev row, prev col)
                    cell.x = row00[i - dx].x + cx;
                    cell.y = row00[i - dx].y + cy;
                    cell.z = row00[i - dx].z + cz;

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
    void DanielssonSlice(int w, int h, CellDelta3 deltas[], int cw)
    {
        CellDelta3* row0 = deltas;
        CellDelta3 cell;

        for (int j = 1; j < h; j++)
        {
            CellDelta3* row1 = row0 + w;

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
            CellDelta3* row1 = row0 - w;

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

void DOL::Danielsson(int w, int h, int d, CellDelta3 deltas[], int cw)
{
    int sliceStride = w * h;

    CellDelta3* slice0 = deltas;
    CellDelta3 cell;

    for (int k = 1; k < d; k++)
    {
        CellDelta3* slice1 = slice0 + sliceStride;

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
        CellDelta3* slice1 = slice0 - sliceStride;

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
    void JumpFlood(int step, int w, int h, int d, const CellDelta3 deltasIn[], CellDelta3 deltasOut[], int cw)
    {
        int wh = w * h;

        #pragma omp parallel for
        for (int k = 0; k < d; k++)
        {
            const CellDelta3* row    = deltasIn  + k * wh;
            CellDelta3*       rowOut = deltasOut + k * wh;

            int kc0 = (k - step) >= 0 ? -step : 0;
            int kc1 = (k + step) <  d ? +step : 0;

            for (int j = 0; j < h; j++)
            {
                int jc0 = (j - step) >= 0 ? -step : 0;
                int jc1 = (j + step) <  h ? +step : 0;

                for (int i = 0; i < w; i++)
                {
                    CellDelta3 minCell = row[i];

                    int ic0 = (i - step) >= 0 ? -step : 0;
                    int ic1 = (i + step) <  w ? +step : 0;

                    for (int kc = kc0; kc <= kc1; kc += step)
                    for (int jc = jc0; jc <= jc1; jc += step)
                    for (int ic = ic0; ic <= ic1; ic += step)
                    {
                        CellDelta3 cell = row[i + ic + jc * w + kc * wh];

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

void DOL::JumpFlood(int w, int h, int d, CellDelta3 deltas[], int cw)
{
    CellDelta3* temp = new CellDelta3[w * h * d];
    CellDelta3* d0 = deltas;
    CellDelta3* d1 = temp;

    int wh = w >= h ? w : h;
    int step = FloorPow2((wh >= d ? wh : d) - 1);

    for ( ; step > 0; step /= 2)
    {
        ::JumpFlood(step, w, h, d, d0, d1, cw);

        CellDelta3* dt = d0;
        d0 = d1;
        d1 = dt;
    }

    if (d0 != deltas)
        ::JumpFlood(1, w, h, d, d0, deltas, cw);    // we could copy, but might as well do another round

    delete[] temp;
}

// Brute force -- for checking only!
void DOL::BruteForce(int w, int h, int d, CellDelta3 deltas[], int cw)
{
    for (int k = 0; k < d; k++)
    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++)
    {
        CellDelta3* cell = deltas + k * w * h + j * w + i;

        if (abs(cell->x) >= cw || abs(cell->y) >= cw || abs(cell->z) >= cw) // non-border cell
            continue;

        // For each occupied cell, iterate over all other cells and check for minimal distance
        for (int z = 0; z < d; z++)
        for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            CellDelta3 d3 = { int16_t((i - x) * cw + cell->x), int16_t((j - y) * cw + cell->y), int16_t((k - z) * cw + cell->z) };
            CellDelta3* cell2 = deltas + z * w * h + y * w + x;

            if (Closer(*cell2, d3))
                *cell2 = d3;
        }
    }
}
