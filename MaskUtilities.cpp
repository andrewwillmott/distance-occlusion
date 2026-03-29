//
// MaskUtilities.cpp
//
// Utilities for dealing with 1-bit images/volumes
//
// Andrew Willmott
//

#include "MaskUtilities.hpp"

#include <string.h>

#if defined(DOS_DEBUG) || defined(DEBUG)
    #include <assert.h>
    #define DOL_ASSERT(e) assert(e)
#else
    #define DOL_ASSERT(e) ((void) 0)
#endif

using namespace DOL;

// Declarations

namespace
{
    inline int RNG(uint32_t& seed, int limit)
    {
        seed = uint32_t(seed * uint64_t(1103515245)) + 12345;
        return int((seed * (uint64_t(limit))) >> 32);
    }

    inline void Swap(float& a, float& b) { float t = a; a = b; b = t; }
}


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

void DOL::ClearBitMask(int w, int h, uint32_t mask[])
{
    memset(mask, 0, MaskSize(w, h) * sizeof(uint32_t));
}

void DOL::BitMaskAddPoints(int w, int h, int count, uint32_t seed, uint32_t mask[])
{
    int ws = MaskSize(w);

    for (int i = 0; i < count; i++)
    {
        int x = RNG(seed, w);
        int y = RNG(seed, h);

        DOL_ASSERT(x < w);
        DOL_ASSERT(y < h);

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
        DOL_ASSERT(n < 32);
        uint32_t bit = (1 << n);
        return ~(bit - 1);
    }

    uint32_t RightMask(int n)
    {
        DOL_ASSERT(n < 32);
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

        DOL_ASSERT(x < w);
        DOL_ASSERT(y < h);
        DOL_ASSERT(z < d);

        int fx = x & 31;
        x >>= 5;

        int index = z * h * ws + y * ws + x;
        mask[index] |= 1 << fx;
    }
}

void DOL::ClearBitMask(int w, int h, int d, uint32_t mask[])
{
    memset(mask, 0, MaskSize(w, h, d) * sizeof(uint32_t));
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


//
// 2D
//

void DOL::CreateBitMaskFromTriangles
(
    int             triCount,
    const int       indices[],
    const Vec2f     vertices[],
    const Bounds2&  bbox,
    int w, int h,
    uint32_t        mask[]
)
{
    MaskDraw2 maskDraw(bbox, w, h, mask);

    for (int iv = 0; iv < triCount * 3; iv += 3)
    {
        int iv1 = iv + 0;
        int iv2 = iv + 1;
        int iv3 = iv + 2;

        if (indices)
        {
            iv1 = indices[iv1];
            iv2 = indices[iv2];
            iv3 = indices[iv3];
        }

        Vec2f p1 = vertices[iv1];
        Vec2f p2 = vertices[iv2];
        Vec2f p3 = vertices[iv3];

        maskDraw.DrawTriangle(p1, p2, p3);
    }
}

MaskDraw2::MaskDraw2(const Bounds2&  bbox, int w_in, int h_in, uint32_t mask_in[]) : mask(mask_in), w(w_in), h(h_in)
{
    rowStride = (w + 31) / 32;

    wh = { float(w), float(h) };

    cellMin  = bbox.mMin;
    cellW    = (bbox.mMax - bbox.mMin) / wh;
    cellInvW = wh / (bbox.mMax - bbox.mMin);

    hw = cellW * 0.500001f;
}

void MaskDraw2::DrawTriangle(Vec2f p1, Vec2f p2, Vec2f p3, float conservative)
{
    Bounds2 triBounds = TriangleBounds(p1, p2, p3);

    Vec2f cellsMin = (triBounds.mMin - cellMin) * cellInvW;
    Vec2f cellsMax = (triBounds.mMax - cellMin) * cellInvW;

    int cx0 = vl_max(FloorToInt32(cellsMin.x), 0);
    int cy0 = vl_max(FloorToInt32(cellsMin.y), 0);

    int cx1 = vl_min( CeilToInt32(cellsMax.x), w);
    int cy1 = vl_min( CeilToInt32(cellsMax.y), h);

    if (cx0 == cx1)
        return;

    if (cx0 + 1 == cx1 && cy0 + 1 == cy1)   // quick out for single-cell triangle
    {
        size_t cellBits = rowStride * cy0 + (cx0 >> 5);
        mask[cellBits] |= 1 << (cx0 & 0x1F);
        return;
    }

    Vec2f centre = {};
    int rowBits = rowStride * cy0 + (cx0 >> 5);

    TriCellIntersectDelta2 triCell(p1, p2, p3, hw * conservative);

    for (int y = cy0; y < cy1; y++)
    {
        int cellBits = rowBits;
        rowBits += rowStride;

        centre.y = cellMin.y + cellW.y * (y + 0.5f);

        if (!triCell.IntersectY(centre))
            continue;

        int i = (cx0 & 0x1F);

        for (int x = cx0; x < cx1; x++)
        {
            DOL_ASSERT(cellBits >= 0 && cellBits < w * h);

            centre.x = cellMin.x + cellW.x * (x + 0.5f);

            if (triCell.IntersectX(centre))
                mask[cellBits] |= 1 << i;

            if (++i == 32)
            {
                cellBits++;
                i = 0;
            }
        }
    }
}

TriCellIntersectDelta2::TriCellIntersectDelta2(Vec2f p1, Vec2f p2, Vec2f p3, Vec2f hwIn) :
    hw(hwIn)
{
    e1 = p2 - p1;
    e2 = p3 - p2;
    e3 = p1 - p3;

    Vec2f vmin, vmax;
    vmin = -hw;
    vmax =  hw;

    Vec2f ae1 = abs(e1);
    Vec2f ae2 = abs(e2);
    Vec2f ae3 = abs(e3);

    t1 =  e1.y * p2.x - e1.x * p2.y;
    u1 =  e1.y * p3.x - e1.x * p3.y;

    if (t1 >= u1) Swap(t1, u1);

    t2 =  e2.y * p1.x - e2.x * p1.y;
    u2 =  e2.y * p2.x - e2.x * p2.y;

    if (t2 >= u2) Swap(t2, u2);

    t3 =  e3.y * p2.x - e3.x * p2.y;
    u3 =  e3.y * p3.x - e3.x * p3.y;

    if (t3 >= u3) Swap(t3, u3);

    float r1, r2, r3;
    r1 = ae1.y * hw.x + ae1.x * hw.y;
    r2 = ae2.y * hw.x + ae2.x * hw.y;
    r3 = ae3.y * hw.x + ae3.x * hw.y;

    t1 -= r1;
    t2 -= r2;
    t3 -= r3;

    u1 += r1;
    u2 += r2;
    u3 += r3;

    pMin = MinElts(p1, p2);
    pMin = MinElts(pMin, p3);

    pMax = MaxElts(p1, p2);
    pMax = MaxElts(pMax, p3);
}

bool TriCellIntersectDelta2::IntersectX(Vec2f c)
{
    float dr;

    dr =  e1.y * c.x - e1.x * c.y;

    if (t1 > dr || u1 < dr)
        return false;

    dr =  e2.y * c.x - e2.x * c.y;

    if (t2 > dr || u2 < dr)
        return false;

    dr =  e3.y * c.x - e3.x * c.y;

    if (t3 > dr || u3 < dr)
        return false;

    if (pMin.x - c.x > hw.x || pMax.x - c.x < -hw.x)
        return false;

    return true;
}

//
// 3D
//

void DOL::CreateBitMaskFromTriangles
(
    int             triCount,
    const int       indices[],
    const Vec3f     vertices[],
    const Bounds3&  bbox,
    int w, int h, int d,
    uint32_t        mask[]
)
{
    MaskDraw3 maskDraw(bbox, w, h, d, mask);

    for (int iv = 0; iv < triCount * 3; iv += 3)
    {
        int iv1 = iv + 0;
        int iv2 = iv + 1;
        int iv3 = iv + 2;

        if (indices)
        {
            iv1 = indices[iv1];
            iv2 = indices[iv2];
            iv3 = indices[iv3];
        }

        Vec3f p1 = vertices[iv1];
        Vec3f p2 = vertices[iv2];
        Vec3f p3 = vertices[iv3];

        maskDraw.DrawTriangle(p1, p2, p3);
    }
}

MaskDraw3::MaskDraw3(const Bounds3& bbox, int w_in, int h_in, int d_in, uint32_t mask_in[]) :
    mask(mask_in), w(w_in), h(h_in), d(d_in)
{
    rowStride    = (w + 31) / 32;
    sliceStride  = rowStride * h;

    whd      = { float(w), float(h), float(d) };

    cellMin  = bbox.mMin;
    cellW    = (bbox.mMax - bbox.mMin) / whd;
    cellInvW = whd / (bbox.mMax - bbox.mMin);

    hw = cellW * 0.500001f;
}

void MaskDraw3::DrawTriangle(Vec3f p1, Vec3f p2, Vec3f p3, float conservative)
{
    Bounds3 triBounds = TriangleBounds(p1, p2, p3);

    Vec3f cellsMin = (triBounds.mMin - cellMin) * cellInvW;
    Vec3f cellsMax = (triBounds.mMax - cellMin) * cellInvW;

    int cx0 = vl_max(FloorToInt32(cellsMin.x), 0);
    int cy0 = vl_max(FloorToInt32(cellsMin.y), 0);
    int cz0 = vl_max(FloorToInt32(cellsMin.z), 0);

    int cx1 = vl_min( CeilToInt32(cellsMax.x), w);
    int cy1 = vl_min( CeilToInt32(cellsMax.y), h);
    int cz1 = vl_min( CeilToInt32(cellsMax.z), d);

    if (cx0 == cx1)
        return;

    if (cx0 + 1 == cx1 && cy0 + 1 == cy1 && cz0 + 1 == cz1)   // quick out for single-cell triangle
    {
        size_t cellBits = sliceStride * cz0 + rowStride * cy0 + (cx0 >> 5);
        mask[cellBits] |= 1 << (cx0 & 0x1F);
        return;
    }

    Vec3f centre = { 0, 0, 0 };
    int sliceBits = sliceStride * cz0 + rowStride * cy0 + (cx0 >> 5);

    TriCellIntersectDelta3 triCell(p1, p2, p3, hw * conservative);

    for (int z = cz0; z < cz1; z++)
    {
        int rowBits = sliceBits;
        sliceBits += sliceStride;

        centre.z = cellMin.z + cellW.z * (z + 0.5f);

        if (!triCell.IntersectZ(centre))
            continue;

        for (int y = cy0; y < cy1; y++)
        {
            int cellBits = rowBits;
            rowBits += rowStride;

            centre.y = cellMin.y + cellW.y * (y + 0.5f);

            if (!triCell.IntersectY(centre))
                continue;

            int i = (cx0 & 0x1F);

            for (int x = cx0; x < cx1; x++)
            {
                DOL_ASSERT(cellBits >= 0 && cellBits < w * h * d);

                centre.x = cellMin.x + cellW.x * (x + 0.5f);

                if (triCell.IntersectX(centre))
                    mask[cellBits] |= 1 << i;

                if (++i == 32)
                {
                    cellBits++;
                    i = 0;
                }
            }
        }
    }
}

TriCellIntersectDelta3::TriCellIntersectDelta3(Vec3f p1, Vec3f p2, Vec3f p3, Vec3f hwIn) :
    hw(hwIn)
{
    e1 = p2 - p1;
    e2 = p3 - p2;
    e3 = p1 - p3;

    normal = TriangleDoubleAreaNormal(p1, p2, p3);

    Vec3f vmin, vmax;

    for (int i = 0; i < 3; i++)
        if ((&normal.x)[i] > 0.0f)
        {
            vmin[i] = -hw[i];
            vmax[i] = +hw[i];
        }
        else
        {
            vmin[i] = +hw[i];
            vmax[i] = -hw[i];
        }

    nDotVMin = dot(normal, vmin);
    nDotVMax = dot(normal, vmax);
    nd       = dot(normal, p1);

    Vec3f ae1 = abs(e1);
    Vec3f ae2 = abs(e2);
    Vec3f ae3 = abs(e3);

    t1.x =  e1.z * p1.y - e1.y * p1.z;
    t1.y = -e1.z * p1.x + e1.x * p1.z;
    t1.z =  e1.y * p2.x - e1.x * p2.y;

    u1.x =  e1.z * p3.y - e1.y * p3.z;
    u1.y = -e1.z * p3.x + e1.x * p3.z;
    u1.z =  e1.y * p3.x - e1.x * p3.y;

    if (t1.x >= u1.x) Swap(t1.x, u1.x);
    if (t1.y >= u1.y) Swap(t1.y, u1.y);
    if (t1.z >= u1.z) Swap(t1.z, u1.z);

    t2.x =  e2.z * p1.y - e2.y * p1.z;
    t2.y = -e2.z * p1.x + e2.x * p1.z;
    t2.z =  e2.y * p1.x - e2.x * p1.y;

    u2.x =  e2.z * p3.y - e2.y * p3.z;
    u2.y = -e2.z * p3.x + e2.x * p3.z;
    u2.z =  e2.y * p2.x - e2.x * p2.y;

    if (t2.x >= u2.x) Swap(t2.x, u2.x);
    if (t2.y >= u2.y) Swap(t2.y, u2.y);
    if (t2.z >= u2.z) Swap(t2.z, u2.z);

    t3.x =  e3.z * p1.y - e3.y * p1.z;
    t3.y = -e3.z * p1.x + e3.x * p1.z;
    t3.z =  e3.y * p2.x - e3.x * p2.y;

    u3.x =  e3.z * p2.y - e3.y * p2.z;
    u3.y = -e3.z * p2.x + e3.x * p2.z;
    u3.z =  e3.y * p3.x - e3.x * p3.y;

    if (t3.x >= u3.x) Swap(t3.x, u3.x);
    if (t3.y >= u3.y) Swap(t3.y, u3.y);
    if (t3.z >= u3.z) Swap(t3.z, u3.z);

    Vec3f r1, r2, r3;
    r1.x = ae1.z * hw.y + ae1.y * hw.z;
    r1.y = ae1.z * hw.x + ae1.x * hw.z;
    r1.z = ae1.y * hw.x + ae1.x * hw.y;

    r2.x = ae2.z * hw.y + ae2.y * hw.z;
    r2.y = ae2.z * hw.x + ae2.x * hw.z;
    r2.z = ae2.y * hw.x + ae2.x * hw.y;

    r3.x = ae3.z * hw.y + ae3.y * hw.z;
    r3.y = ae3.z * hw.x + ae3.x * hw.z;
    r3.z = ae3.y * hw.x + ae3.x * hw.y;

    t1 -= r1;
    t2 -= r2;
    t3 -= r3;

    u1 += r1;
    u2 += r2;
    u3 += r3;

    pMin = MinElts(p1, p2);
    pMin = MinElts(pMin, p3);

    pMax = MaxElts(p1, p2);
    pMax = MaxElts(pMax, p3);
}

bool TriCellIntersectDelta3::IntersectY(const Vec3f& c)
{
    float dr;

    dr =  e1.z * c.y - e1.y * c.z;

    if (t1.x > dr || u1.x < dr)
        return false;

    dr =  e2.z * c.y - e2.y * c.z;

    if (t2.x > dr || u2.x < dr)
        return false;

    dr =  e3.z * c.y - e3.y * c.z;

    if (t3.x > dr || u3.x < dr)
        return false;

    if (pMin.y - c.y > hw.y || pMax.y - c.y < -hw.y)
        return false;

    return true;
}

bool TriCellIntersectDelta3::IntersectX(const Vec3f& c)
{
    float dr;

    dr = -e1.z * c.x + e1.x * c.z;

    if (t1.y > dr || u1.y < dr)
        return false;

    dr =  e1.y * c.x - e1.x * c.y;

    if (t1.z > dr || u1.z < dr)
        return false;

    dr = -e2.z * c.x + e2.x * c.z;

    if (t2.y > dr || u2.y < dr)
        return false;

    dr =  e2.y * c.x - e2.x * c.y;

    if (t2.z > dr || u2.z < dr)
        return false;

    dr = -e3.z * c.x + e3.x * c.z;

    if (t3.y > dr || u3.y < dr)
        return false;

    dr =  e3.y * c.x - e3.x * c.y;

    if (t3.z > dr || u3.z < dr)
        return false;

    if (pMin.x - c.x > hw.x || pMax.x - c.x < -hw.x)
        return false;

    // Finally test against plane. Timing tests confirm it's best to have this last.
    float d = nd - dot(normal, c);
    return (nDotVMin <= d) && (nDotVMax >= d);
}
