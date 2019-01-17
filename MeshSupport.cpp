//------------------------------------------------------------------------------
// Purpose: Quick and dirty mesh support
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#define  _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include "MeshSupport.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef _MSC_VER
    #pragma warning (disable: 4244)
    #define strcasecmp _stricmp
    #define strtok_r strtok_s
#endif

using namespace MSL;

// Finds size of volume bbox given allowed error 'eps'. Derivation from Malmer et al.
Bounds3f MSL::FindAOBounds(float eps, Bounds3f modelBounds)
{
    float s = 1.0f / (4.0f * float(M_PI) * eps);

    Vec3f d =
    {
        modelBounds.mMax.x - modelBounds.mMin.x,
        modelBounds.mMax.y - modelBounds.mMin.y,
        modelBounds.mMax.z - modelBounds.mMin.z
    };

    Vec3f a =
    {
        sqrtf(d.y * d.z * s),
        sqrtf(d.z * d.x * s),
        sqrtf(d.x * d.y * s)
    };

    Bounds3f result =
    {
        {
            modelBounds.mMin.x - a.x,
            modelBounds.mMin.y - a.y,
            modelBounds.mMin.z - a.z,
        },
        {
            modelBounds.mMax.x + a.x,
            modelBounds.mMax.y + a.y,
            modelBounds.mMax.z + a.z,
        }
    };

    return result;
}

// This version allows you to supply your own axial projected area bounds
Bounds3f MSL::FindAOBounds(float eps, Bounds3f modelBounds, Vec3f minArea, Vec3f maxArea)
{
    float s = 1.0f / (4.0f * float(M_PI) * eps);

    Bounds3f result =
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

namespace MSL
{
    inline Vec3f  operator- (Vec3f  a, Vec3f b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
    inline Vec3f  operator+ (Vec3f  a, Vec3f b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
    inline Vec3f  operator* (Vec3f  a, Vec3f b) { return { a.x * b.x, a.y * b.y, a.z * b.z }; }
    inline Vec3f  operator/ (Vec3f  a, Vec3f b) { return { a.x / b.x, a.y / b.y, a.z / b.z }; }
    inline Vec3f& operator+=(Vec3f& a, Vec3f b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
    inline Vec3f& operator-=(Vec3f& a, Vec3f b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }

    inline float dot(Vec3f a, Vec3f b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    inline Vec3f abs(Vec3f v)          { return { fabsf(v.x), fabsf(v.y), fabsf(v.z) }; }
    inline float len(Vec3f v)          { return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z); }

    template <typename T> inline void swap(T& a, T& b) { T t(a); a = b; b = t; }

    template<class T> inline T Max(T a, T b)
    {
        return b < a ? a : b;
    }

    template<class T> inline T Min(T a, T b)
    {
        return a < b ? a : b;
    }

    inline Vec3f MinElts(const Vec3f& a, const Vec3f& b)
    {
        return Vec3f
        {
            Min(a.x, b.x),
            Min(a.y, b.y),
            Min(a.z, b.z)
        };
    }

    inline Vec3f MaxElts(const Vec3f& a, const Vec3f& b)
    {
        return Vec3f
        {
            Max(a.x, b.x),
            Max(a.y, b.y),
            Max(a.z, b.z)
        };
    }

    inline int FloorToInt32(float x)
    {
        return int(floorf(x));
    }
    inline int CeilToInt32(float x)
    {
        return int(ceilf(x));
    }

    Vec3f TriDoubleAreaNormal
    (
        const Vec3f& a,
        const Vec3f& b,
        const Vec3f& c
    )
    {
        Vec3f n =
        {
            (a.y - b.y) * (a.z + b.z) + (b.y - c.y) * (b.z + c.z) + (c.y - a.y) * (c.z + a.z),
            (a.z - b.z) * (a.x + b.x) + (b.z - c.z) * (b.x + c.x) + (c.z - a.z) * (c.x + a.x),
            (a.x - b.x) * (a.y + b.y) + (b.x - c.x) * (b.y + c.y) + (c.x - a.x) * (c.y + a.y)
        };

        return n;
    }

    struct cTriCellIntersectDelta
    /// Helper class for intersecting a triangle with many identically-sized,
    /// axis-aligned cells, as in voxelisation. It is assumed the cells will
    /// be visited in z/y/x order, and helpers are provided for early outs where
    /// an entire slice or row is guaranteed not to intersect.
    {
        cTriCellIntersectDelta(Vec3f p1, Vec3f p2, Vec3f p3, Vec3f hw); ///< Initialise from triangle and cell half-width

        bool IntersectZ(const Vec3f& c);    ///< Returns true if cell extended infinitely in z would intersect (slice).
        bool IntersectY(const Vec3f& c);    ///< Returns true if cell extended infinitely in y would intersect (row), assuming IntersectDeltaZ().
        bool IntersectX(const Vec3f& c);    ///< Returns true if cell would intersect, assuming IntersectDeltaZ() && IntersectDeltaY().

        // Pre-calculated triangle intersection data
        Vec3f e1, t1, u1;   // Triangle edges + axis testing info for axis testing
        Vec3f e2, t2, u2;
        Vec3f e3, t3, u3;

        Vec3f pMin;         // For range testing
        Vec3f pMax;
        Vec3f hw;           // half cell width

        Vec3f normal;       // For plane testing
        float nDotVMin;
        float nDotVMax;
        float nd;
    };

    cTriCellIntersectDelta::cTriCellIntersectDelta(Vec3f p1, Vec3f p2, Vec3f p3, Vec3f hwIn) :
        hw(hwIn)
    {
        e1 = p2 - p1;
        e2 = p3 - p2;
        e3 = p1 - p3;

        normal = TriDoubleAreaNormal(p1, p2, p3);

        Vec3f vmin, vmax;

        for (int i = 0; i < 3; i++)
            if ((&normal.x)[i] > 0.0f)
            {
                (&vmin.x)[i] = -(&hw.x)[i];
                (&vmax.x)[i] = +(&hw.x)[i];
            }
            else
            {
                (&vmin.x)[i] = +(&hw.x)[i];
                (&vmax.x)[i] = -(&hw.x)[i];
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

        if (t1.x >= u1.x) swap(t1.x, u1.x);
        if (t1.y >= u1.y) swap(t1.y, u1.y);
        if (t1.z >= u1.z) swap(t1.z, u1.z);

        t2.x =  e2.z * p1.y - e2.y * p1.z;
        t2.y = -e2.z * p1.x + e2.x * p1.z;
        t2.z =  e2.y * p1.x - e2.x * p1.y;

        u2.x =  e2.z * p3.y - e2.y * p3.z;
        u2.y = -e2.z * p3.x + e2.x * p3.z;
        u2.z =  e2.y * p2.x - e2.x * p2.y;

        if (t2.x >= u2.x) swap(t2.x, u2.x);
        if (t2.y >= u2.y) swap(t2.y, u2.y);
        if (t2.z >= u2.z) swap(t2.z, u2.z);

        t3.x =  e3.z * p1.y - e3.y * p1.z;
        t3.y = -e3.z * p1.x + e3.x * p1.z;
        t3.z =  e3.y * p2.x - e3.x * p2.y;

        u3.x =  e3.z * p2.y - e3.y * p2.z;
        u3.y = -e3.z * p2.x + e3.x * p2.z;
        u3.z =  e3.y * p3.x - e3.x * p3.y;

        if (t3.x >= u3.x) swap(t3.x, u3.x);
        if (t3.y >= u3.y) swap(t3.y, u3.y);
        if (t3.z >= u3.z) swap(t3.z, u3.z);

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

    inline bool cTriCellIntersectDelta::IntersectZ(const Vec3f& c)
    {
        if (pMin.z - c.z > hw.z || pMax.z - c.z < -hw.z)
            return false;

        return true;
    }

    bool cTriCellIntersectDelta::IntersectY(const Vec3f& c)
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

    bool cTriCellIntersectDelta::IntersectX(const Vec3f& c)
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

    inline Bounds3f TriangleBounds(const Vec3f& a, const Vec3f& b, const Vec3f& c)
    {
        Bounds3f box;

        box.mMin = MinElts(a, b);
        box.mMin = MinElts(box.mMin, c);

        box.mMax = MaxElts(a, b);
        box.mMax = MaxElts(box.mMax, c);

        return box;
    }
}

void MSL::CreateBitMaskFromTriangles
(
    int             triCount,
    const int       indices[],
    const Vec3f     vertices[],
    const Bounds3f& bbox,
    int w, int h, int d,
    uint32_t        mask[]
)
{
    int rowStride    = (w + 31) / 32;
    int sliceStride  = rowStride * h;

    Vec3f whd      = { float(w), float(h), float(d) };

    Vec3f cellMin  = bbox.mMin;
    Vec3f cellW    = (bbox.mMax - bbox.mMin) / whd;
    Vec3f cellInvW = whd / (bbox.mMax - bbox.mMin);

    Vec3f hw = cellW * Vec3f{0.500001f, 0.500001f, 0.500001f};

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

        Vec3f p1 = vertices[indices[iv + 0]];
        Vec3f p2 = vertices[indices[iv + 1]];
        Vec3f p3 = vertices[indices[iv + 2]];

        Bounds3f triBounds  = TriangleBounds(p1, p2, p3);

        Vec3f cellsMin = (triBounds.mMin - cellMin) * cellInvW;
        Vec3f cellsMax = (triBounds.mMax - cellMin) * cellInvW;

        int cx0 = Max(FloorToInt32(cellsMin.x), 0);
        int cy0 = Max(FloorToInt32(cellsMin.y), 0);
        int cz0 = Max(FloorToInt32(cellsMin.z), 0);

        int cx1 = Min( CeilToInt32(cellsMax.x), w);
        int cy1 = Min( CeilToInt32(cellsMax.y), h);
        int cz1 = Min( CeilToInt32(cellsMax.z), d);

        if (cx0 == cx1)
            continue;

        if (cx0 + 1 == cx1 && cy0 + 1 == cy1 && cz0 + 1 == cz1)   // quick out for single-cell triangle
        {
            size_t cellBits = sliceStride * cz0 + rowStride * cy0 + (cx0 >> 5);
            mask[cellBits] |= 1 << (cx0 & 0x1F);
            continue;
        }

        Vec3f centre = { 0, 0, 0 };
        int sliceBits = sliceStride * cz0 + rowStride * cy0 + (cx0 >> 5);

        cTriCellIntersectDelta triCell(p1, p2, p3, hw);

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
                    assert(cellBits >= 0 && cellBits < w * h * d);

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
}

void MSL::CreateDirW8FromTriangles
(
    int             triCount,
    const int       indices[],
    const Vec3f     vertices[],
    const Bounds3f& bbox,
    int w, int h, int d,
    float           dirW8[]
)
{
    int rowStride    = w;
    int sliceStride  = rowStride * h;
    int volumeStride = sliceStride * d;

    Vec3f whd      = { float(w), float(h), float(d) };

    Vec3f cellMin  = bbox.mMin;
    Vec3f cellW    = (bbox.mMax - bbox.mMin) / whd;
    Vec3f cellInvW = whd / (bbox.mMax - bbox.mMin);

    Vec3f hw = cellW * Vec3f{0.500001f, 0.500001f, 0.500001f};

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

        Vec3f p1 = vertices[indices[iv + 0]];
        Vec3f p2 = vertices[indices[iv + 1]];
        Vec3f p3 = vertices[indices[iv + 2]];

        Bounds3f triBounds  = TriangleBounds(p1, p2, p3);

        Vec3f cellsMin = (triBounds.mMin - cellMin) * cellInvW;
        Vec3f cellsMax = (triBounds.mMax - cellMin) * cellInvW;

        int cx0 = Max(FloorToInt32(cellsMin.x), 0);
        int cy0 = Max(FloorToInt32(cellsMin.y), 0);
        int cz0 = Max(FloorToInt32(cellsMin.z), 0);

        int cx1 = Min( CeilToInt32(cellsMax.x), w);
        int cy1 = Min( CeilToInt32(cellsMax.y), h);
        int cz1 = Min( CeilToInt32(cellsMax.z), d);

        if (cx0 == cx1)
            continue;

        cTriCellIntersectDelta triCell(p1, p2, p3, hw);

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
            float bc =  (i & 1) ? -triCell.normal.x : triCell.normal.x;
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

                            assert(dirIndex >= 0 && dirIndex < volumeStride * 8);

                            if (dirW8[dirIndex] < dirValue[i])
                                dirW8[dirIndex] = dirValue[i];
                        }
                    }
                }
            }
        }
    }
}


namespace
{
    void InPlaceTriangulate(int numVerts, std::vector<int>& indices)
    {
        // Assume polygon of numVerts indices at the end of 'indices', and triangulate it in place.

        int baseVertexIndex = (int) indices.size() - numVerts;
        int baseEltIndex = indices[baseVertexIndex];

        for (int i = 1; i < numVerts - 2; i++)
        {
            indices.insert(indices.begin() + baseVertexIndex + 3 * i, indices[baseVertexIndex + 3 * i - 1]);
            indices.insert(indices.begin() + baseVertexIndex + 3 * i, baseEltIndex);
        }
    }

    bool FaceCommand(cMesh* mesh, int argc, const char* va[])
    {
        argc--;
        va++;

        char* end;

        for (int i = 0; i < argc; i++)
        {
            long ip = strtol(va[i], &end, 10);

            if (end != va[i])
                mesh->mPositionIndices.push_back(int(ip) - 1);

            if (end[0] == '/')
            {
                const char* next = end + 1;
                long it = strtol(next, &end, 10);

                if (end != next)
                    mesh->mUVIndices.push_back(int(it) - 1);
            }

            if (end[0] == '/')
            {
                const char* next = end + 1;
                long in = strtol(end + 1, &end, 10);

                if (end != next)
                    mesh->mNormalIndices.push_back(int(in) - 1);
            }
        }

        // quick and dirty in-place triangulation.
        if (argc > 3)
        {
            InPlaceTriangulate(argc, mesh->mPositionIndices);

            if (!mesh->mNormalIndices.empty())
                InPlaceTriangulate(argc, mesh->mNormalIndices);

            if (!mesh->mUVIndices.empty())
                InPlaceTriangulate(argc, mesh->mUVIndices);
        }

        return true;
    }

    bool PositionCommand(cMesh* mesh, int argc, const char* va[])
    {
        if (argc < 4)
            return false;

        Vec3f p;
        p.x = atof(va[1]);
        p.y = atof(va[2]);
        p.z = atof(va[3]);

        mesh->mPositions.push_back(p);

        return true;
    }

    bool NormalCommand(cMesh* mesh, int argc, const char* va[])
    {
        if (argc < 4)
            return false;

        Vec3f p;
        p.x = atof(va[1]);
        p.y = atof(va[2]);
        p.z = atof(va[3]);

        mesh->mNormals.push_back(p);

        return true;
    }

    bool TexCoordCommand(cMesh* mesh, int argc, const char* va[])
    {
        if (argc < 3)
            return false;

        Vec2f p;
        p.x = atof(va[1]);
        p.y = atof(va[2]);

        mesh->mUVs.push_back(p);

        return true;
    }

    bool ObjectCommand(cMesh*, int, const char**)
    {
        return true;
    }

    bool GroupCommand(cMesh*, int, const char**)
    {
        return true;
    }

    bool SmoothingGroupCommand(cMesh*, int, const char**)
    {
        return true;
    }

    bool MaterialCommand(cMesh*, int, const char**)
    {
        return true;
    }
    bool MaterialLibraryCommand(cMesh*, int, const char**)
    {
        return true;
    }


    bool ProcessObjCommand(cMesh* mesh, int argc, const char* argv[])
    {
        assert(argc > 0);

        switch (argv[0][0])
        {
            case 'f':
                return FaceCommand(mesh, argc, argv);
            case 'v':
                if (argv[0][1] == 'n')
                    return NormalCommand(mesh, argc, argv);
                else if (argv[0][1] == 't')
                    return TexCoordCommand(mesh, argc, argv);
                else if (argv[0][1] == 0)
                    return PositionCommand(mesh, argc, argv);
                break;
            case 'o':
                return ObjectCommand(mesh, argc, argv);
            case 'g':
                return GroupCommand(mesh, argc, argv);
            case 's':
                return SmoothingGroupCommand(mesh, argc, argv);
            case 'u':
                if (strcasecmp(argv[0], "usemtl") == 0)
                    return MaterialCommand(mesh, argc, argv);
                break;
            case 'm':
                if (strcasecmp(argv[0], "mtllib") == 0)
                    return MaterialLibraryCommand(mesh, argc, argv);
                break;
            case '#':
                return true;
        }
        
        return false;
    }


    int Split(char* buffer, int maxArgs, const char** argv)
    {
        char* last = 0;
        maxArgs--;  // always reserve the last spot for 0 terminator
        
        for (int argc = 0; argc < maxArgs; argc++)
        {
            argv[argc] = strtok_r(buffer, " \t\n\r", &last);
            buffer = 0;

            if (!argv[argc] || !argv[argc][0])
                return argc;
        }

        fprintf(stderr, "Warning: ignored arguments past %d\n", maxArgs - 1);
        argv[maxArgs] = 0;
        return maxArgs;
    }
}

bool MSL::ReadObjFile(FILE* file, cMesh* mesh)
{
    *mesh = cMesh();

    const int kMaxArgs = 256;
    const char* argv[kMaxArgs];

    char lineBuffer[1024];

    while (fgets(lineBuffer, 1024, file))
    {
        int argc = Split(lineBuffer, kMaxArgs, argv);

        if (argc > 0)
        {
            if (!ProcessObjCommand(mesh, argc, argv))
            {
                fprintf(stderr, "Can't parse command: '%s'\n", argv[0]);
                return false;
            }
        }
    }

    printf("Read %zu positions, %zu normals, %zu uvs\n",
        mesh->mPositions.size(),
        mesh->mNormals.size(),
        mesh->mUVs.size()
    );

    printf("     %zd position indices, %zd normal indices, %zd uv indices\n",
        mesh->mPositionIndices.size(),
        mesh->mNormalIndices.size(),
        mesh->mUVIndices.size()
    );

    return true;
}

Bounds3f MSL::FindBounds(const cMesh& mesh)
{
    Bounds3f bbox = { { +FLT_MAX, +FLT_MAX, +FLT_MAX }, { -FLT_MAX, -FLT_MAX, -FLT_MAX } };

    for (int i = 0, n = (int) mesh.mPositions.size(); i < n; i++)
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

    return bbox;
}


