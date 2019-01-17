//------------------------------------------------------------------------------
// Purpose: Distance Field and Occlusion Generation
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#ifndef DISTANCE_OCCLUSION_H
#define DISTANCE_OCCLUSION_H

#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

namespace DOL
{
    // --- Bit Masks -----------------------------------------------------------

    // We store input occlusion as 1-bit bitmasks, padded to word boundaries on
    // each row. E.g., a 100 x 100 image would be 100 x 4 uint32_ts.

    int MaskSize(int w);                 // return mask size in words
    int MaskSize(int w, int h);          // return mask size in words
    int MaskSize(int w, int h, int d);   // return mask size in words

    uint32_t* CreateBitMask (int w, int h);         // Create an empty image bit mask
    uint32_t* CreateBitMask (int w, int h, int d);  // Create an empty volume bit mask
    void      DestroyBitMask(uint32_t*& mask);      // Free the given mask

    // Routines for generating simple 1-bit 'occlusion' test images/volumes
    void BitMaskAddPoints(int w, int h,        int n, uint32_t seed, uint32_t mask[]);
    void BitMaskAddPoints(int w, int h, int d, int n, uint32_t seed, uint32_t mask[]);
    // Add 'n' random points to the bit mask

    void BitMaskAddBlock(int w, int h,        int sides, uint32_t mask[]);
    void BitMaskAddBlock(int w, int h, int d, int sides, uint32_t mask[]);
    // Add a box with the given number of sides to the bit mask


    // --- Cell deltas ---------------------------------------------------------

    struct cCellDelta2
    {
        int16_t x;
        int16_t y;
    };
    // Stores delta from current cell to closest solid boundary

    struct cCellDelta3
    {
        int16_t x;
        int16_t y;
        int16_t z;
    };
    // Stores delta from current cell to closest solid boundary

    const int16_t kMaxDelta2 = INT16_MAX;
    const int16_t kMaxDelta3 = 25000;       // INT16_MAX will overflow 32-bit arithmetic when calculating distances, reduce max to fix.

    // --- Distance ------------------------------------------------------------

    void InitDistancesFromBitMask(int w, int h, const uint32_t mask[], int16_t distances[], int16_t maxD = INT16_MAX);
    void InitDistancesFromBitMask(int w, int h, const uint32_t mask[], int32_t distances[], int32_t maxD = INT32_MAX);
    void InitDistancesFromBitMask(int w, int h, const uint32_t mask[], float   distances[], float   maxD = FLT_MAX);
    // Given a mask, initialises a distance field to 0 where bits are one, INT16_MAX/INT32_MAX otherwise

    int Chamfer(int w, int h, int16_t distances[]);
    int Chamfer(int w, int h, int32_t distances[]);
    // Use Chamfer algorithm to propagate distances directly. Not as accurate as the delta-based approaches below.
    // Returns multiplier -- divide through by this to get distance in cells.

    void Felzenszwalb(int w, int h, float distances[]);
    // Use Felzenszwalb/Huttenlocher algorithm to propagate squared distances via parabolas.

    // Delta-based algorithms
    void InitDeltasFromBitMask(int w, int h,        const uint32_t mask[], cCellDelta2 delta[], bool invert = false);
    void InitDeltasFromBitMask(int w, int h, int d, const uint32_t mask[], cCellDelta3 delta[], bool invert = false);
    // Given a mask, initialises a distance field to 0 where bits are one, kMaxDelta otherwise

    void InitBorderDeltasFromBitMask(int w, int h, const uint32_t mask[], cCellDelta2 delta[]);
    void InitBorderDeltasFromBitMask(int w, int h, int d, const uint32_t mask[], cCellDelta3 delta[]);
    // Given a mask, initialises a distance field correctly at 0/1 boundaries, fills the remainder with kMaxDelta.
    // This is the variant to use to generate signed distance fields (in conjunction with the original mask to indicate sign.)

    void FastSweep (int w, int h,        cCellDelta2 delta[], int cw = 1);
    void FastSweep (int w, int h, int d, cCellDelta3 delta[], int cw = 1);
    // Find final distance field from initial distance field using fast sweep algorithm.

    void Danielsson(int w, int h,        cCellDelta2 delta[], int cw = 1);
    void Danielsson(int w, int h, int d, cCellDelta3 delta[], int cw = 1);
    // Use Danielsson algorithm to propagate cell deltas -- similar to FastSweep but combines passes.

    void JumpFlood (int w, int h,        cCellDelta2 delta[], int cw = 1);
    void JumpFlood (int w, int h, int d, cCellDelta3 delta[], int cw = 1);
    // Use "jump flooding" -- intended as a parallel algorithm, here for result comparison.

    void BruteForce(int w, int h,        cCellDelta2 delta[], int cw = 1);
    void BruteForce(int w, int h, int d, cCellDelta3 delta[], int cw = 1);
    // Finds optimal delta values using brute force. Useful for checking results of faster algorithms.

    // --- Occlusion -----------------------------------------------------------

    void InitDirWFromBitMask(int w, int h,        const uint32_t mask[], float dirW4[]);
    void InitDirWFromBitMask(int w, int h, int d, const uint32_t mask[], float dirW8[]);
    // Initialise 4-way or 8-way direction field from given bitmask.

    void OcclusionSweep(int w, int h,        float dirW4[]);
    // Calculate occlusion in each quadrant's direction. dirW4 should have been initialised with some InitDirW variant.
    void OcclusionSweep(int w, int h, int d, float dirW8[]);
    void OcclusionSweep(int w, int h, int d, const float inDirW8[], float outDirW8[]);
    // Stores occlusion caused by inDirW in outDirW in each octant's direction. (Hence can be used to avoid self occlusion.)

    uint8_t EncodeAO4  (const float dirW[4]);   // Encode 2D direction components into overall occlusion.
    uint8_t EncodeAO8  (const float dirW[8]);   // Encode 3D direction components into overall occlusion.
    void    EncodeDirW4(const float dirW[4], uint8_t aoDir[4]); // Encode 2D direction components into directed occlusion in x/y and overall occlusion in z.
    void    EncodeDirW8(const float dirW[8], uint8_t aoDir[4]); // Encode 3D direction components into directed occlusion in x/y/z and overall occlusion in a.


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

    inline uint8_t EncodeU8(float v)
    {
        v = (v < 0.0f) ? 0.0f : ((v > 1.0f) ? 1.0f : v);
        return (uint8_t) lrintf(v * 255);
    }

    inline uint8_t EncodeAO4(const float dirW[4])
    {
        float aoF = dirW[0] + dirW[1] + dirW[2] + dirW[3];
        
        return EncodeU8(1.0f - 0.25f * aoF);
    }

    inline uint8_t EncodeAO8(const float dirW[8])
    {
        float aoF = dirW[0] + dirW[1] + dirW[2] + dirW[3] + dirW[4] + dirW[5] + dirW[6] + dirW[7];
        
        return EncodeU8(1.0f - 0.125f * aoF);
    }

    inline void EncodeDirW4(const float dirW[4], uint8_t aoDir[4])
    {
        float aoFW = dirW[0] + dirW[1] + dirW[2] + dirW[3];
        float aoFX = dirW[0] - dirW[1] + dirW[2] - dirW[3];
        float aoFY = dirW[0] + dirW[1] - dirW[2] - dirW[3];
        
        float invLen = 1.0f / sqrtf(aoFX * aoFX + aoFY * aoFY);

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

    inline void EncodeDirW8(const float dirW[8], uint8_t aoDir[4])
    {
        float aoFW = dirW[0] + dirW[1] + dirW[2] + dirW[3] + dirW[4] + dirW[5] + dirW[6] + dirW[7];
        float aoFX = dirW[0] - dirW[1] + dirW[2] - dirW[3] + dirW[4] - dirW[5] + dirW[6] - dirW[7];
        float aoFY = dirW[0] + dirW[1] - dirW[2] - dirW[3] + dirW[4] + dirW[5] - dirW[6] - dirW[7];
        float aoFZ = dirW[0] + dirW[1] + dirW[2] + dirW[3] - dirW[4] - dirW[5] - dirW[6] - dirW[7];
        
        float invLen = 1.0f / sqrtf(aoFX * aoFX + aoFY * aoFY + aoFZ * aoFZ);
        
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
}

#endif
