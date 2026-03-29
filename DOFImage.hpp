//
// DOFImage.hpp
//
// Distance/occlusion field image IO
//
// Andrew Willmott
//

#ifndef DOF_IMAGE_H
#define DOF_IMAGE_H

#include <stdint.h>

namespace DOL
{
    // Image helpers for debugging

    uint32_t* ReadMask (const char* path, int threshold, bool lessThan, int& w, int& h);
    void      WriteMask(const char* name, int w, int h, const uint32_t mask[]);

    // Distances

    void WriteDistances(const char* name, const char* type, int w, int h, float maxLen, const float[]);
    void WriteDistances(const char* name, const char* type, int w, int h, float maxLen, const int16_t[]);

    struct CellDelta2;
    struct CellDelta3;

    void WriteDistances(const char* name, const char* type, int w, int h, float maxLen, const CellDelta2 deltas[]);
    void WriteDistances(const char* name, const char* type, int w, int h, int d, float maxLen, const CellDelta3 deltas[]);

    // Signed distances (mask + delta or delta in, delta out)

    enum DFImageType : int
    {
        kImageStandard,
        kImageDeltas,
        kImageRG,
        kImageAngles,
    };

    void WriteDistances(const char* name, int w, int h, float maxLen, const uint32_t   mask   [], const CellDelta2 deltas [], int imageType);
    void WriteDistances(const char* name, int w, int h, float maxLen, const CellDelta2 deltas0[], const CellDelta2 deltas1[], int imageType);

    void WriteDistances(const char* name, int w, int h, int d, float maxLen, const uint32_t   mask   [], const CellDelta3 deltas [], int imageType);
    void WriteDistances(const char* name, int w, int h, int d, float maxLen, const CellDelta3 deltas0[], const CellDelta3 deltas1[], int imageType);

    // Occlusion
    void WriteOcclusion(const char* name, int w, int h,        const float dirW[]);
    void WriteOcclusion(const char *name, int w, int h, int d, const float dirW[]);

    void WriteShadow(const char* name, const char* type, int w, int h, const float shadow[], int index = -1);
}

#endif
