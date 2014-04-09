//------------------------------------------------------------------------------
// Purpose: Distance Field and Occlusion Generation
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#ifndef DISTANCE_OCCLUSION_H
#define DISTANCE_OCCLUSION_H

#include <stdint.h>
#include <stdlib.h>

namespace DOL
{
    struct cCellDelta2
    {
        int16_t x;
        int16_t y;
    };
    ///< Stores delta from current cell to closest solid cell

    struct cCellDelta3
    {
        int16_t x;
        int16_t y;
        int16_t z;
    };
    ///< Stores delta from current cell to closest solid cell

    // --- Test ----------------------------------------------------------------

    // Routines for generating simple 1-bit 'occlusion' test images/volumes
    void CreateBitMask(int w, int h,        uint32_t seed, uint32_t mask[]);
    void CreateBitMask(int w, int h, int d, uint32_t seed, uint32_t mask[]);
    ///< Creates a random occlusion mask with the given dimensions

    void CreateBitMaskFromBlock(int w, int h,        int sides, uint32_t mask[]);
    void CreateBitMaskFromBlock(int w, int h, int d, int sides, uint32_t mask[]);
    ///< Create a bitmask from a box with the given number of sides

    int MaskSize(int w);                 // return mask size in words
    int MaskSize(int w, int h);          // return mask size in words
    int MaskSize(int w, int h, int d);   // return mask size in words

    // --- Distance ------------------------------------------------------------

    void InitDFFromBitMask(int w, int h, const uint32_t mask[], float distances[]);
    void InitDFFromBitMask(int w, int h, const uint32_t mask[], int32_t distances[]);
    void InitDeltaFromBitMask(int w, int h, const uint32_t mask[], cCellDelta2 delta[]);
    void InitDeltaFromBitMask(int w, int h, int d, const uint32_t mask[], cCellDelta3 delta[]);
    ///< Given a mask, initialises a distance field to 0 where bits are one, FLT_MAX/INT32_MAX otherwise

    void Chamfer(int w, int h, float distances[]);
    ///< Use Chamfer algorithm to propagate float distances directly
    void Chamfer(int w, int h, int32_t distances[]);
    ///< Use Chamfer algorithm to propagate int distances directly

    void Danielsson(int w, int h,        cCellDelta2 delta[]);
    void Danielsson(int w, int h, int d, cCellDelta3 delta[]);
    ///< Use Danielsson algorithm to propagate cell deltas.

    void FastSweep(int w, int h,        cCellDelta2 delta[]);
    void FastSweep(int w, int h, int d, cCellDelta3 delta[]);
    ///< Find final distance field from initial distance field using fast sweep algorithm.

    void FindOptimalDeltas(int w, int h,        cCellDelta2 delta[]);
    void FindOptimalDeltas(int w, int h, int d, cCellDelta3 delta[]);
    ///< Finds optimal delta values using brute force. Useful for checking results of faster algorithms.

    // --- Occlusion -----------------------------------------------------------

    void InitDirWFromBitMask(int w, int h,        const uint32_t mask[], float dirW4[]);
    void InitDirWFromBitMask(int w, int h, int d, const uint32_t mask[], float dirW8[]);
    ///< Initialise 4-way direction field from given bitmask.

    void OcclusionSweep(int w, int h,        float dirW4[]);
    ///< Calculate occlusion in each quadrant's direction. dirW4 should have been initialised with some InitDirW variant.
    void OcclusionSweep(int w, int h, int d, float dirW8[]);
    void OcclusionSweep(int w, int h, int d, const float inDirW8[], float outDirW8[]);
    ///< Stores occlusion caused by inDirW in outDirW in each octant's direction. (Hence can be used to avoid self occlusion.)


    // --- Inlines -------------------------------------------------------------

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
}

#endif
