//
// DistanceFields.hpp
//
// Distance field generation utilities
//
// Andrew Willmott
//

#ifndef DISTANCE_FIELDS_H
#define DISTANCE_FIELDS_H

#include "MaskUtilities.hpp"

#include <float.h>

namespace DOL
{
    // --- Cell deltas ---------------------------------------------------------

    struct CellDelta2
    {
        int16_t x;
        int16_t y;
    };
    // Stores delta from current cell to closest solid boundary

    struct CellDelta3
    {
        int16_t x;
        int16_t y;
        int16_t z;
    };
    // Stores delta from current cell to closest solid boundary

    const int16_t kMaxDelta2 = INT16_MAX;
    const int16_t kMaxDelta3 = 25000;       // INT16_MAX will overflow 32-bit arithmetic when calculating distances, reduce max to fix.

    int32_t SqrDistance(const CellDelta2& cd);
    int32_t SqrDistance(const CellDelta3& cd);
    float      Distance(const CellDelta2& cd);
    float      Distance(const CellDelta3& cd);

    // --- Distance ------------------------------------------------------------

    void InitDistancesFromBitMask(int w, int h, const uint32_t mask[], int16_t distances[], int16_t maxD = INT16_MAX);
    void InitDistancesFromBitMask(int w, int h, const uint32_t mask[], int32_t distances[], int32_t maxD = INT32_MAX);
    void InitDistancesFromBitMask(int w, int h, const uint32_t mask[], float   distances[], float   maxD = FLT_MAX);
    // Given a mask, initialises a distance field to 0 where bits are one, INT16_MAX/INT32_MAX otherwise

    void InitDistancesFromBitMask(int w, int h, const uint32_t mask[], float   distances[], float minD, float maxD);  // set minD where bits are 1, maxD otherwise.

    int Chamfer(int w, int h, int16_t distances[]);
    int Chamfer(int w, int h, int32_t distances[]);
    // Use Chamfer algorithm to propagate distances directly. Not as accurate as the delta-based approaches below.
    // Returns multiplier -- divide through by this to get distance in cells.

    void Felzenszwalb(int w, int h, float distances[]);
    // Use Felzenszwalb/Huttenlocher algorithm to propagate squared distances via parabolas.

    // Delta-based algorithms
    void InitDeltasFromBitMask(int w, int h,        const uint32_t mask[], CellDelta2 delta[], bool invert = false);
    void InitDeltasFromBitMask(int w, int h, int d, const uint32_t mask[], CellDelta3 delta[], bool invert = false);
    // Given a mask, initialises a distance field to 0 where bits are one, kMaxDelta otherwise

    void InitBorderDeltasFromBitMask(int w, int h, const uint32_t mask[], CellDelta2 delta[]);
    void InitBorderDeltasFromBitMask(int w, int h, int d, const uint32_t mask[], CellDelta3 delta[]);
    // Given a mask, initialises a distance field correctly at 0/1 boundaries, fills the remainder with kMaxDelta.
    // This is the variant to use to generate signed distance fields (in conjunction with the original mask to indicate sign.)

    void FastSweep (int w, int h,        CellDelta2 delta[], int cw = 1);
    void FastSweep (int w, int h, int d, CellDelta3 delta[], int cw = 1);
    // Find final distance field from initial distance field using fast sweep algorithm.

    void Danielsson(int w, int h,        CellDelta2 delta[], int cw = 1);
    void Danielsson(int w, int h, int d, CellDelta3 delta[], int cw = 1);
    // Use Danielsson algorithm to propagate cell deltas -- similar to FastSweep but combines passes.

    void JumpFlood (int w, int h,        CellDelta2 delta[], int cw = 1);
    void JumpFlood (int w, int h, int d, CellDelta3 delta[], int cw = 1);
    // Use "jump flooding" -- intended as a parallel algorithm, here for result comparison.

    void BruteForce(int w, int h,        CellDelta2 delta[], int cw = 1);
    void BruteForce(int w, int h, int d, CellDelta3 delta[], int cw = 1);
    // Finds optimal delta values using brute force. Useful for checking results of faster algorithms.


    // --- Inlines -------------------------------------------------------------

    inline int32_t SqrDistance(const CellDelta2& cd)
    {
        return cd.x * cd.x + cd.y * cd.y;
    }

    inline int32_t SqrDistance(const CellDelta3& cd)
    {
        return cd.x * cd.x + cd.y * cd.y + cd.z * cd.z;
    }

    inline float Distance(const CellDelta2& cd)
    {
        return sqrtf((float) SqrDistance(cd));
    }

    inline float Distance(const CellDelta3& cd)
    {
        return sqrtf((float) SqrDistance(cd));
    }
}

#endif
