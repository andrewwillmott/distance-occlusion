//
// OcclusionFields.cpp
//
// Andrew Willmott
//

#include "OcclusionFields.hpp"

#include "MaskUtilities.hpp"

#if defined(DOS_DEBUG) || defined(DEBUG)
    #include <assert.h>
    #define DOL_ASSERT assert
#else
    #define DOL_ASSERT(e) ((void) 0)
#endif

using namespace DOL;

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
                    DOL_ASSERT(32 * j + i < w);

                    block[i + 0 * s] = 1.0f;
                    block[i + 1 * s] = 1.0f;
                    block[i + 2 * s] = 1.0f;
                    block[i + 3 * s] = 1.0f;
                }

                m >>= 1;
            }
        }
}

void DOL::InitDirWFromBitMask(int w, int h, int d, const uint32_t mask[], float dirW[])
// Initialise octants from bit mask.
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
    inline float QuadrantCombine(float a, float b)
    {
        // Find new quadrant estimate given same quad estimates from neighbour cells in sweep direction.

        constexpr float kQC1 = 0.46f;   // value with one neighbour fully occluded
        constexpr float kQC2 = 0.975f;  // value with both neighbours fully occluded

        return kQC1 * (a + b) + (kQC2 - 2 * kQC1) * a * b;
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

namespace
{
    inline float OctantCombine(float a, float b, float c)
    {
        // Find new octant estimate given same octant estimates from neighbour cells in sweep direction.

        constexpr float kOC1 = 0.3333333f;  // value with one neighbour fully occluded
        constexpr float kOC3 = 0.6666667f;  // value with all neighbours fully occluded

        // return kOC1 * (a + b + c) + (kOC3 - 3 * kOC1) * (a * b + b * c + c * a - a * b * c) * 0.5f;

        // Cheaper, works okay with smaller kOC3.
        return kOC1 * (a + b + c) + (kOC3 - 3 * kOC1) * (a * b * c);
    }

    // ix = input occlusion, ox = output occlusion
    inline float OctantCombine(float oa, float ob, float oc, float ia, float ib, float ic)
    {
        float a = ia + (1.0f - ia) * oa;
        float b = ib + (1.0f - ib) * ob;
        float c = ic + (1.0f - ic) * oc;

        return OctantCombine(a, b, c);
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
    int stride = w * h * d;

    for (int sweep = 0; sweep < 8; sweep++)
    {
        const int sx = (sweep >> 0) & 1;
        const int sy = (sweep >> 1) & 1;
        const int sz = (sweep >> 2) & 1;

        ::OcclusionSweep(w, h, d, sx, sy, sz, dirW + sweep * stride);
    }
}

void DOL::OcclusionSweep(int w, int h, int d, const float* occDirW, float* dirW)
{
    int stride = w * h * d;

    for (int sweep = 0; sweep < 8; sweep++)
    {
        const int sx = (sweep >> 0) & 1;
        const int sy = (sweep >> 1) & 1;
        const int sz = (sweep >> 2) & 1;

        ::OcclusionSweep(w, h, d, sx, sy, sz, occDirW + sweep * stride, dirW + sweep * stride);
    }
}

void DOL::FillDirAOFromDirW4(int w, int h, const float dirW4[], RGBA8 aoDir[])
{
    int sliceStride = w * h;

    const float* sDirW[4];
    for (int i = 0; i < 4; i++)
        sDirW[i] = dirW4 + sliceStride * i;

    for (int i = 0, n = sliceStride; i < n; i++)
    {
        float aoFW = sDirW[0][i] + sDirW[1][i] + sDirW[2][i] + sDirW[3][i];
        aoFW *= 0.25f;

        float aoFX = sDirW[0][i] - sDirW[1][i] + sDirW[2][i] - sDirW[3][i];
        float aoFY = sDirW[0][i] + sDirW[1][i] - sDirW[2][i] - sDirW[3][i];

        float invLen = 1.0f / (sqrtf(aoFX * aoFX + aoFY * aoFY) + 1e-8f);

        // sum of all solid angles
        aoFX = aoFW * invLen * aoFX * 0.5f + 0.5f;
        aoFY = aoFW * invLen * aoFY * 0.5f + 0.5f;

        aoDir->c[0] = RoundToInt32(aoFX * 255);
        aoDir->c[1] = RoundToInt32(aoFY * 255);
        aoDir->c[2] = 0;
        aoDir->c[3] = RoundToInt32(aoFW * 255);

        aoDir++;
    }
}

void DOL::FillDirAOFromDirW8(int w, int h, int d, const float dirW8[], RGBA8 aoDir[])
{
    int sliceStride = w * h;
    int volStride   = sliceStride * d;

    for (int slice = 0; slice < d; slice++)
    {
        int sliceOffset = sliceStride * slice;

        const float* sDirW[8];
        for (int i = 0; i < 8; i++)
            sDirW[i] = dirW8 + sliceOffset + volStride * i;

        for (int i = 0, n = sliceStride; i < n; i++)
        {
            float aoFW = sDirW[0][i] + sDirW[1][i] + sDirW[2][i] + sDirW[3][i]
                       + sDirW[4][i] + sDirW[5][i] + sDirW[6][i] + sDirW[7][i];
            aoFW *= 0.125;

            float aoFX = sDirW[0][i] - sDirW[1][i] + sDirW[2][i] - sDirW[3][i]
                       + sDirW[4][i] - sDirW[5][i] + sDirW[6][i] - sDirW[7][i];
            float aoFY = sDirW[0][i] + sDirW[1][i] - sDirW[2][i] - sDirW[3][i]
                       + sDirW[4][i] + sDirW[5][i] - sDirW[6][i] - sDirW[7][i];
            float aoFZ = sDirW[0][i] + sDirW[1][i] + sDirW[2][i] + sDirW[3][i]
                       - sDirW[4][i] - sDirW[5][i] - sDirW[6][i] - sDirW[7][i];

            float invLen = 1.0f / (sqrtf(aoFX * aoFX + aoFY * aoFY + aoFZ * aoFZ) + 1e-8f);

            // sum of all solid angles
            aoFX = aoFW * invLen * aoFX * 0.5f + 0.5f;
            aoFY = aoFW * invLen * aoFY * 0.5f + 0.5f;
            aoFZ = aoFW * invLen * aoFZ * 0.5f + 0.5f;

            aoDir->c[0] = RoundToInt32(aoFX * 255);
            aoDir->c[1] = RoundToInt32(aoFY * 255);
            aoDir->c[2] = RoundToInt32(aoFZ * 255);
            aoDir->c[3] = RoundToInt32(aoFW * 255);

            aoDir++;
        }
    }
}

// Finds size of volume bbox given allowed error 'eps'. Derivation from Malmer et al.
Bounds3 DOL::FindAOBounds(float eps, Bounds3 modelBounds)
{
    float s = 1.0f / (4.0f * vlf_pi * eps);

    Vec3f d = modelBounds.mMax - modelBounds.mMin;

    Vec3f a =
    {
        sqrtf(d.y * d.z * s),
        sqrtf(d.z * d.x * s),
        sqrtf(d.x * d.y * s)
    };

    Bounds3 result =
    {
        modelBounds.mMin - a,
        modelBounds.mMax + a
    };

    return result;
}

// This version allows you to supply your own axial projected area bounds
Bounds3 DOL::FindAOBounds(float eps, Bounds3 modelBounds, Vec3f minArea, Vec3f maxArea)
{
    float s = 1.0f / (4.0f * vlf_pi * eps);

    Bounds3 result =
    {
        {
            modelBounds.mMin.x - sqrtf(minArea.x * s),
            modelBounds.mMin.y - sqrtf(minArea.y * s),
            modelBounds.mMin.z - sqrtf(minArea.z * s)
        },
        {
            modelBounds.mMax.x + sqrtf(maxArea.x * s),
            modelBounds.mMax.y + sqrtf(maxArea.y * s),
            modelBounds.mMax.z + sqrtf(maxArea.z * s)
        }
    };

    return result;
}

void DOL::CreateDirW8FromTriangles
(
    int            triCount,
    const int      indices[],
    const Vec3f    vertices[],
    const Bounds3& bbox,
    int w, int h, int d,
    float          dirW8[]
)
{
    int rowStride    = w;
    int sliceStride  = h * rowStride;
    int volumeStride = d * sliceStride;

    Vec3f whd      = { float(w), float(h), float(d) };

    Vec3f cellMin  = bbox.mMin;
    Vec3f cellW    = (bbox.mMax - bbox.mMin) / whd;
    Vec3f cellInvW = whd / (bbox.mMax - bbox.mMin);

    Vec3f hw = cellW * 0.500001f;

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
            continue;

        TriCellIntersectDelta3 triCell(p1, p2, p3, hw);

        float area = 0.5f * len(triCell.normal);

        // We have two conflicting desires here -- we want the reconstruction of the occlusion direction
        // to be the same as the triangle normal, and ideally the occlusion amount = 0.5 = the hemisphere.
        // Luckily for our basis M, MtM = 8I, and more so, if you zero negative weights,
        // MtM = 4I. As we treat each weight as an octant coverage value however, we have an issue
        // in that e_z -> b_i = 1, and e_111 -> b_i = (3, 1, 1, 1) / sqrt(3), whose max coverage > 1, which
        // our sweep algorithm doesn't handle.

        float bs = 0.577350269f / area; //  1 / ||a|| sqrt(3)

        int   dirOffset[8];
        float dirValue [8];
        int   dirCount = 0;

        for (int i = 0; i < 8; i++)
        {
            float bc;
            bc =  (i & 1) ? -triCell.normal.x : triCell.normal.x;
            bc += (i & 2) ? -triCell.normal.y : triCell.normal.y;
            bc += (i & 4) ? -triCell.normal.z : triCell.normal.z;

            if (bc > 1e-6f)
            {
                dirOffset[dirCount] = i * volumeStride;
                dirValue [dirCount] = bc * bs;
                dirCount++;
            }
        }

        Vec3f centre = { 0, 0, 0 };

        for (int z = cz0; z < cz1; z++)
        {
            centre.z = cellMin.z + cellW.z * (z + 0.5f);

            if (!triCell.IntersectZ(centre))
                continue;

            for (int y = cy0; y < cy1; y++)
            {
                centre.y = cellMin.y + cellW.y * (y + 0.5f);

                if (!triCell.IntersectY(centre))
                    continue;

                for (int x = cx0; x < cx1; x++)
                {
                    centre.x = cellMin.x + cellW.x * (x + 0.5f);

                    if (triCell.IntersectX(centre))
                    {
                        int index = z * sliceStride + y * rowStride + x;

                        for (int i = 0; i < dirCount; i++)
                        {
                            int dirIndex = dirOffset[i] + index;

                            DOL_ASSERT(dirIndex >= 0 && dirIndex < volumeStride * 8);

                            if (dirW8[dirIndex] < dirValue[i])
                                dirW8[dirIndex] = dirValue[i];
                        }
                    }
                }
            }
        }
    }
}
float DOL::BlobRadiusFromEpsilon(float eps)
{
    float x = sqrtf(vlf_pi * eps);
    return x / (1.0f + x);
}

namespace
{
    inline float BFromSphere(float d, float r)    // normalised solid angle (0-2) of sphere of radius r from distance d
    {
        r /= d;

        if (r > 1.0f)
            return 1.0f;

        return 1 - sqrtf(1 - r * r);
    }

    inline Vec4f DirBFromSphere(Vec3f v, float r)
    {
        float id = 1.0f / len(v);
        float rd = id * r;
        float b = (rd > 1.0f) ? 1 : 1 - sqrtf(1 - rd * rd);

        Vec3f dirB = (b * id) * v;
        return { dirB.x, dirB.y, dirB.z, b };
    }

#ifdef UNUSED
    inline float WFromSphere(float d, float r)
    {
        return vlf_twoPi * BFromSphere(d, r);
    }
    inline Vec4f DirWFromSphere(Vec3f v, float r)
    {
        return vlf_twoPi * DirBFromSphere(v, r);
    }
#endif
}

void DOL::FillDirAOWithBlob(int w, int h, int d, float r, uint8_t* aoDir)
{
    Vec3f s(2.0f / w, 2.0f / h, 2.0f / d);

    float sMax = vl_max(vl_max(s.x, s.y), s.z);
    float maxEdge = BFromSphere(1.0f - 0.5f * sMax, r);
    float rescale = 1.0f / (1.0f - 2.0f * maxEdge);

    for (int z = 0; z < d; z++)
    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        Vec4f dirB = DirBFromSphere(Vec3f((x + 0.5f) * s.x - 1.0f, (y + 0.5f) * s.y - 1.0f, (z + 0.5f) * s.z - 1.0f), r);

        dirB.w = (dirB.w - maxEdge) * rescale;
        if (dirB.w < 0.0f)
            dirB.w = 0.0f;

        float aoFW = dirB.w * 0.5f;
        float aoFX = dirB.x * 0.25f + 0.5f;
        float aoFY = dirB.y * 0.25f + 0.5f;
        float aoFZ = dirB.z * 0.25f + 0.5f;

        aoDir[0] = RoundToInt32(aoFX * 255.0f);
        aoDir[1] = RoundToInt32(aoFY * 255.0f);
        aoDir[2] = RoundToInt32(aoFZ * 255.0f);
        aoDir[3] = RoundToInt32(aoFW * 255.0f);

        aoDir += 4;
    }
}

void DOL::TaperVolume(int w, int h, int d, RGBA8* dirW)
{
    int rowStride   = w;
    int sliceStride = w * h;
    int volStride   = w * h * d;

    uint8_t maxEdgeW = 0;

    const RGBA8* planeW = dirW;            // z min plane
    for (int i = 0; i < sliceStride; i++)
        if (maxEdgeW < planeW[i].a)
            maxEdgeW = planeW[i].a;

    planeW = dirW + volStride - sliceStride;    // z max plane
    for (int i = 0; i < sliceStride; i++)
        if (maxEdgeW < planeW[i].a)
            maxEdgeW = planeW[i].a;

    planeW = dirW;                              // y min plane
    for (int z = 1; z < d - 1; z++)
    {
        const RGBA8* rowW = planeW + z * sliceStride;

        for (int x = 0; x < w; x++)
            if (maxEdgeW < rowW[x].a)
                maxEdgeW = rowW[x].a;
    }

    planeW = dirW + sliceStride - rowStride;    // y max plane
    for (int z = 1; z < d - 1; z++)
    {
        const RGBA8* rowW = planeW + z * sliceStride;

        for (int x = 0; x < w; x++)
            if (maxEdgeW < rowW[x].a)
                maxEdgeW = rowW[x].a;
    }

    planeW = dirW + rowStride;                  // x min/max planes
    for (int z = 1; z < d - 1; z++)
    {
        const RGBA8* rowW0 = planeW + z * sliceStride;
        const RGBA8* rowW1 = rowW0 + rowStride - 1;

        for (int y = 1; y < h - 1; y++)
        {
            if (maxEdgeW < rowW0->a)
                maxEdgeW = rowW0->a;

            if (maxEdgeW < rowW1->a)
                maxEdgeW = rowW1->a;

            rowW0 += rowStride;
            rowW1 += rowStride;
        }
    }

    if (maxEdgeW > 128)
        maxEdgeW = 128;

    uint8_t t = maxEdgeW;
    float s = 255.0f / (255.0f - t);

    for (int i = 0; i < volStride; i++)
    {
        float sw = s * (dirW[i].a - t);
        if (sw < 0.0f)
            sw = 0.0f;

        DOL_ASSERT(sw <= 255.0f);
        dirW[i].a = RoundToInt32(sw);
    }
}

void DOL::EncodeDirW4(const float dirW[4], uint8_t aoDir[4])
{
    float aoFW = dirW[0] + dirW[1] + dirW[2] + dirW[3];
    float aoFX = dirW[0] - dirW[1] + dirW[2] - dirW[3];
    float aoFY = dirW[0] + dirW[1] - dirW[2] - dirW[3];

    float invLen = 1.0f / (sqrtf(aoFX * aoFX + aoFY * aoFY) + 1e-8f);

    aoFX *= invLen;
    aoFY *= invLen;
    aoFW *= 0.25f;

    aoFX = aoFW * aoFX * 0.5f + 0.5f;
    aoFY = aoFW * aoFY * 0.5f + 0.5f;

    aoDir[0] = EncodeU8(aoFX);
    aoDir[1] = EncodeU8(aoFY);
    aoDir[2] = EncodeU8(aoFW);
    aoDir[3] = 255;
}

void DOL::EncodeDirW8(const float dirW[8], uint8_t aoDir[4])
{
    float aoFW = dirW[0] + dirW[1] + dirW[2] + dirW[3] + dirW[4] + dirW[5] + dirW[6] + dirW[7];
    float aoFX = dirW[0] - dirW[1] + dirW[2] - dirW[3] + dirW[4] - dirW[5] + dirW[6] - dirW[7];
    float aoFY = dirW[0] + dirW[1] - dirW[2] - dirW[3] + dirW[4] + dirW[5] - dirW[6] - dirW[7];
    float aoFZ = dirW[0] + dirW[1] + dirW[2] + dirW[3] - dirW[4] - dirW[5] - dirW[6] - dirW[7];

    float invLen = 1.0f / (sqrtf(aoFX * aoFX + aoFY * aoFY + aoFZ * aoFZ) + 1e-8f);

    aoFX *= invLen;
    aoFY *= invLen;
    aoFZ *= invLen;
    aoFW *= 0.125f;

    aoFX = aoFW * aoFX * 0.5f + 0.5f;
    aoFY = aoFW * aoFY * 0.5f + 0.5f;
    aoFZ = aoFW * aoFZ * 0.5f + 0.5f;

    aoDir[0] = EncodeU8(aoFX);
    aoDir[1] = EncodeU8(aoFY);
    aoDir[2] = EncodeU8(aoFZ);
    aoDir[3] = EncodeU8(aoFW);
}


namespace
{
    void ShadowSweepRows(float g, int sx, int sy, int w, int h, const float* occlusion, float* shadow)
    {
        const float g1 = g;
        const float g0 = 1.0f - g1;

        const int dx = 1 - 2 * sx;
        const int dy = 1 - 2 * sy;

        const int di = dx;
        const int ib = sx * (w - 1);
        const int ie = ib + di * w;

        const int dj = dy * w;
        const int jb = sy * (h - 1) * w;

        // Propagation: d' = m + d * (1 - m), where m = lerp of previous row cells

        const float* iRow0 = occlusion + jb;
        float*       oRow0 = shadow    + jb;

        // do first row -- only depends on itself
        for (int i = ib; i != ie; i += di)
            oRow0[i] = iRow0[i];

        for (int j = 1; j < h; j++)
        {
            const float* iRow1 = iRow0 + dj;
            float*       oRow1 = oRow0 + dj;

            // do first cell, only depends on prev row above
            float m = oRow0[ib] * g0;
            oRow1[ib] = m + (1.0f - m) * iRow1[ib];

            // rest of row -- prev cell above and above left
            for (int i = ib + di; i != ie; i += di)
            {
                float m = oRow0[i - di] * g1 + oRow0[i] * g0;
                oRow1[i] = m + (1.0f - m) * iRow1[i];
            }

            iRow0 = iRow1;
            oRow0 = oRow1;
        }
    }

    void ShadowSweepCols(float g, int sx, int sy, int w, int h, const float* occlusion, float* shadow)
    {
        float g1 = g;
        float g0 = 1.0f - g1;

        const int dx = 1 - 2 * sx;
        const int dy = 1 - 2 * sy;

        const int wh = w * h;   // slice stride

        const int ib = (wh - w) * sy;
        const int ie = ib + wh * dy;
        const int di = dy * w;

        const int jb = (w - 1) * sx;
        const int dj = dx;

        const float* iRow0 = occlusion + jb;
        float*       oRow0 = shadow    + jb;

        // do first col -- only depends on itself
        for (int i = ib; i != ie; i += di)
            oRow0[i] = iRow0[i];

        for (int j = 1; j < w; j++)
        {
            const float* iRow1 = iRow0 + dj;
            float*       oRow1 = oRow0 + dj;

            // do first cell, only depends on prev col above
            float m = oRow0[ib] * g0;
            oRow1[ib] = m + (1.0f - m) * iRow1[ib];

            // rest of col -- prev cell above and above left/right
            for (int i = ib + di; i != ie; i += di)
            {
                float m = oRow0[i - di] * g1 + oRow0[i] * g0;
                oRow1[i] = m + (1.0f - m) * iRow1[i];
            }

            iRow0 = iRow1;
            oRow0 = oRow1;
        }
    }
}

void DOL::ShadowSweep(float gx, float gy, int w, int h, const float* occlusion, float* shadow)
{
    const int sx = gx > 0.0f ? 0 : 1;
    const int sy = gy > 0.0f ? 0 : 1;

    if (fabsf(gx) <= fabsf(gy))
        ShadowSweepRows(fabsf(gx / gy), sx, sy, w, h, occlusion, shadow);
    else
        ShadowSweepCols(fabsf(gy / gx), sx, sy, w, h, occlusion, shadow);
}
