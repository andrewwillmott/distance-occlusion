//
// OcclusionFields.hpp
//
// Occlusion field generation utilities
//
// Andrew Willmott
//

#ifndef OCCLUSION_FIELDS_H
#define OCCLUSION_FIELDS_H

#include "VLMini.hpp"
#include <float.h>
#include <stdint.h>

namespace DOL
{
    void InitDirWFromBitMask(int w, int h,        const uint32_t mask[], float dirW4[]);
    void InitDirWFromBitMask(int w, int h, int d, const uint32_t mask[], float dirW8[]);
    // Initialise 4-way or 8-way direction field from given bitmask.

    void OcclusionSweep(int w, int h,        float dirW4[]);  // Calculate occlusion in each quadrant's direction. dirW4 should have been initialised with some InitDirW variant.
    void OcclusionSweep(int w, int h, int d, float dirW8[]);  // Ditto but for octants
    void OcclusionSweep(int w, int h, int d, const float inDirW8[], float outDirW8[]);  // Stores occlusion caused by inDirW in outDirW in each octant's direction. (Hence can be used to avoid self occlusion.)

    void FillDirAOFromDirW4(int w, int h,        const float dirW4[], RGBA8 occDir[]);  // Converts occlusion over quadrants to encoded 2D direction/occlusion
    void FillDirAOFromDirW8(int w, int h, int d, const float dirW8[], RGBA8 occDir[]);  // Converts occlusion over octants to encoded 3D direction/occlusion

    void CreateDirW8FromTriangles
    (
        int             triCount,
        const int       indices[],
        const Vec3f     vertices[],
        const Bounds3&  bbox,
        int w, int h, int d,
        float           dirW8[]
    );
    // Initialises directional occlusion volumes from the given triangle vertices/indices. Indices can be null.
    void TaperVolume(int w, int h, int d, RGBA8 dirOcc[]); // Taper off occlusion at volume edges to avoid discontinuities.

    Bounds3 FindAOBounds(float eps, Bounds3 modelBounds); // Finds size of volume bbox given allowed error 'eps'. Derivation from Malmer et al.
    Bounds3 FindAOBounds(float eps, Bounds3 modelBounds, Vec3f minArea, Vec3f maxArea);  // Variant of FindAOBounds that lets you supply your own axial projected area bounds

    // Analytic versions for testing...
    void  FillDirAOWithBlob(int w, int h, int d, float r, uint8_t aoDir[]);  // Fill aoDir with field from a sphere/ellipsoid of radius r * (w, h, d)
    float BlobRadiusFromEpsilon(float eps);  // Given error epsilon (of the same type as FindAOBounds) find an appropriate blob radius.

    uint8_t EncodeAO4  (const float dirW[4]);   // Encode 2D direction components into overall occlusion.
    uint8_t EncodeAO8  (const float dirW[8]);   // Encode 3D direction components into overall occlusion.
    void    EncodeDirW4(const float dirW[4], uint8_t aoDir[4]); // Encode 2D direction components into directed occlusion in x/y and overall occlusion in z.
    void    EncodeDirW8(const float dirW[8], uint8_t aoDir[4]); // Encode 3D direction components into directed occlusion in x/y/z and overall occlusion in a.

    // 'Swept' directional shadow propagation.
    void ShadowSweep(float gx, float gy, int w, int h, const float occlusion[], float shadow[]);  // Generate shadow (propagated straight-light occlusion)

    template<class SDF, class V> float ConeCast(const SDF& sdf, V o, V d, float w, float tMin = 0.0f, float tMax = FLT_MAX);
    // Simple signed distance field cone cast utility. sdf(p) is assumed to return the distance from p to the nearest surface.
    // Returns the visibility along the ray, where 1.0f is unoccluded, and 0.0f is fully occluded. w is the cone width at unit distance.


    //
    // Inlines
    //

    uint8_t EncodeU8(float v);  // from MaskUtilities.hpp

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

    template<class SDF, class VEC> inline float ConeCast(const SDF& sdf, VEC o, VEC d, float w, float tMin, float tMax)
    {
        float len_d = len(d);
        float invLenD = 1.0f / (len_d + 1e-8f);
        float minStepVis = 1.0f;

        for (float t = tMin; t < tMax; )
        {
            float sd = sdf(o + d * t);
            float cw = w * t * len_d;

            sd += 0.5f * cw;

            if (sd < 1e-2f)
                return 0.0f;

            float stepVis = sd / cw;  // = 1 as surface just touches cone, down to 0 for full hit

            if (minStepVis > stepVis)
                minStepVis = stepVis;

            t += sd * invLenD;
        }

        return minStepVis;
    }
}

#endif
