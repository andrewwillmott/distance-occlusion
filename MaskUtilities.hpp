//
// MaskUtilities.hpp
//
// Utilities for dealing with 1-bit images/volumes
//
// Andrew Willmott
//

#ifndef MASK_UTILITIES_H
#define MASK_UTILITIES_H

#include "VLMini.hpp"
#include <stdint.h>

namespace DOL
{
    struct Bounds2;
    struct Bounds3;

    // --- Bit Masks -----------------------------------------------------------

    // We store input occlusion as 1-bit bitmasks, padded to word boundaries on
    // each row. E.g., a 100 x 100 image would be 100 x 4 uint32_ts.

    typedef uint32_t* BitMask;
    typedef const uint32_t* ConstBitMask;

    int MaskSize(int w);                 // return mask size in words
    int MaskSize(int w, int h);          // return mask size in words
    int MaskSize(int w, int h, int d);   // return mask size in words

    BitMask   CreateBitMask (int w, int h);         // Create an empty image bit mask
    BitMask   CreateBitMask (int w, int h, int d);  // Create an empty volume bit mask
    void      DestroyBitMask(BitMask& mask);        // Free the given mask

    void      ClearBitMask (int w, int h, BitMask mask);         // Create an empty image bit mask
    void      ClearBitMask (int w, int h, int d, BitMask mask);  // Create an empty volume bit mask

    bool MaskBit   (const BitMask mask, int w, int h,        int cx, int cy);
    void SetMaskBit(      BitMask mask, int w, int h,        int cx, int cy, bool bit);
    bool MaskBit   (const BitMask mask, int w, int h, int d, int cx, int cy, int cz);
    void SetMaskBit(      BitMask mask, int w, int h, int d, int cx, int cy, int cz, bool bit);

    // Routines for generating simple 1-bit 'occlusion' test images/volumes
    void BitMaskAddPoints(int w, int h,        int n, uint32_t seed, BitMask mask);
    void BitMaskAddPoints(int w, int h, int d, int n, uint32_t seed, BitMask mask);
    // Add 'n' random points to the bit mask

    void BitMaskAddBlock(int w, int h,        int sides, BitMask mask);
    void BitMaskAddBlock(int w, int h, int d, int sides, BitMask mask);
    // Add a box with the given number of sides to the bit mask

    uint8_t EncodeU8(float v);
    uint8_t EncodeU8Signed(float v);

    // Mesh rasterisation

    void CreateBitMaskFromTriangles
    (
        int            triCount,
        const int      indices[],
        const Vec2f    vertices[],
        const Bounds2& bbox,
        int w, int h,
        uint32_t       mask[]
    );
    // Initialises occupancy mask from the given mesh. 'indices' may be null.

    void CreateBitMaskFromTriangles
    (
        int            triCount,
        const int      indices[],
        const Vec3f    vertices[],
        const Bounds3& bbox,
        int w, int h, int d,
        BitMask        mask
    );
    // Initialises occupancy mask from the given mesh. 'indices' may be null.

    //
    // Helpers for rasterising triangles directly into a mask field.
    //

    struct MaskDraw2
    {
        MaskDraw2(const Bounds2& bbox, int w, int h, BitMask mask);

        void DrawTriangle(Vec2f p1, Vec2f p2, Vec2f p3, float conservative = 1.0f);
        // Rasterise given triangle into 'mask'. If 'conservative' is specified, it selects between
        //  conservative rasterisation (1.0f), cell centres only (0.0f), or a mix thereof.

        uint32_t* mask;
        int w;
        int h;

        int rowStride;
        Vec2f wh;

        Vec2f cellMin;
        Vec2f cellW;
        Vec2f cellInvW;
        Vec2f hw;
    };

    struct MaskDraw3
    {
        MaskDraw3(const Bounds3& bbox, int w, int h, int d, BitMask mask);

        void DrawTriangle(Vec3f p1, Vec3f p2, Vec3f p3, float conservative = 1.0f);
        // Rasterise given triangle into 'mask'. If 'conservative' is specified, it selects between
        //  conservative rasterisation (1.0f), cell centres only (0.0f), or a mix thereof.

        uint32_t* mask;
        int w;
        int h;
        int d;

        int rowStride;
        int sliceStride;

        Vec3f whd;

        Vec3f cellMin;
        Vec3f cellW;
        Vec3f cellInvW;

        Vec3f hw;
    };

    //
    // Underlying rasterisation helpers
    //

    struct TriCellIntersectDelta2
    {
        // Helper class for intersecting a triangle with many identically-sized,
        // axis-aligned cells, as in voxelisation. It is assumed the cells will
        // be visited in z/y/x order, and helpers are provided for early outs where
        // an entire slice or row is guaranteed not to intersect.

        TriCellIntersectDelta2(Vec2f p1, Vec2f p2, Vec2f p3, Vec2f hw); // Initialise from triangle and cell half-width

        bool IntersectY(Vec2f c);    // Returns true if cell extended infinitely in y would intersect (row)
        bool IntersectX(Vec2f c);    // Returns true if cell would intersect given cell, assuming IntersectDeltaY().

        // Pre-calculated triangle intersection data
        Vec2f e1;
        float t1, u1;   // Triangle edges + axis testing info for axis testing
        Vec2f e2;
        float t2, u2;
        Vec2f e3;
        float t3, u3;

        Vec2f pMin;         // For range testing
        Vec2f pMax;
        Vec2f hw;           // half cell width
    };

    struct TriCellIntersectDelta3
    {
        // Helper class for intersecting a triangle with many identically-sized,
        // axis-aligned cells, as in voxelisation. It is assumed the cells will
        // be visited in z/y/x order, and helpers are provided for early outs where
        // an entire slice or row is guaranteed not to intersect.

        TriCellIntersectDelta3(Vec3f p1, Vec3f p2, Vec3f p3, Vec3f hw); // Initialise from triangle and cell half-width

        bool IntersectZ(const Vec3f& c);    // Returns true if cell extended infinitely in z would intersect (slice).
        bool IntersectY(const Vec3f& c);    // Returns true if cell extended infinitely in y would intersect (row), assuming IntersectDeltaZ().
        bool IntersectX(const Vec3f& c);    // Returns true if cell would intersect, assuming IntersectDeltaZ() && IntersectDeltaY().

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

    //
    // Inlines
    //

    inline int MaskSize(int w)
    {
        return (w + 31) / 32;
    }

    inline int MaskSize(int w, int h)
    {
        return ((w + 31) / 32) * h;
    }

    inline int MaskSize(int w, int h, int d)
    {
        return ((w + 31) / 32) * h * d;
    }

    inline bool MaskBit(const BitMask mask, int w, int h, int cx, int cy)
    {
        return (mask[MaskSize(w) * cy + (cx >> 5)] & (1 << (cx & 0x1F))) != 0;
    }

    inline bool MaskBit(const BitMask mask, int w, int h, int d, int cx, int cy, int cz)
    {
        int strideY = MaskSize(w);
        int strideZ = strideY * h;

        return (mask[strideZ * cz + strideY * cy + (cx >> 5)] & (1 << (cx & 0x1F))) != 0;
    }

    inline void SetMaskBit(BitMask mask, int w, int h, int cx, int cy, bool bit)
    {
        if (bit)
            mask[MaskSize(w) * cy + (cx >> 5)] |=  (1 << (cx & 0x1F));
        else
            mask[MaskSize(w) * cy + (cx >> 5)] &= ~(1 << (cx & 0x1F));
    }

    inline void SetMaskBit(BitMask mask, int w, int h, int d, int cx, int cy, int cz, bool bit)
    {
        int strideY = MaskSize(w);
        int strideZ = strideY * h;

        if (bit)
            mask[strideZ * cz + strideY * cy + (cx >> 5)] |=  (1 << (cx & 0x1F));
        else
            mask[strideZ * cz + strideY * cy + (cx >> 5)] &= ~(1 << (cx & 0x1F));
    }

    inline bool TriCellIntersectDelta2::IntersectY(Vec2f c)
    {
        if (pMin.y - c.y > hw.y || pMax.y - c.y < -hw.y)
            return false;

        return true;
    }

    inline bool TriCellIntersectDelta3::IntersectZ(const Vec3f& c)
    {
        if (pMin.z - c.z > hw.z || pMax.z - c.z < -hw.z)
            return false;

        return true;
    }

    inline uint8_t EncodeU8(float v)
    {
        v = (v < 0.0f) ? 0.0f : ((v > 1.0f) ? 1.0f : v);
        return (uint8_t) lrintf(v * 255);
    }

    inline uint8_t EncodeU8Signed(float v)
    {
        v = (v < -1.0f) ? -1.0f : ((v > 1.0f) ? 1.0f : v);
        return (uint8_t) lrintf(128 + v * 127); // ensure 0 = 128
    }
}

#endif
