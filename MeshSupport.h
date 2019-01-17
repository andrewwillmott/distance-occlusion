//------------------------------------------------------------------------------
// Purpose: Quick and dirty mesh support
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#ifndef MESH_SUPPORT_H
#define MESH_SUPPORT_H

#include <vector>
#include <stdio.h>

namespace MSL
{
    // Mesh
    struct Vec2f
    {
        float x;
        float y;
    };

    struct Vec3f
    {
        float x;
        float y;
        float z;
    };

    struct Bounds3f
    {
        Vec3f mMin;
        Vec3f mMax;
    };

    Bounds3f FindBounds();

    Bounds3f FindAOBounds(float eps, Bounds3f modelBounds);
    // Finds size of volume bbox given allowed error 'eps'. Derivation from Malmer et al.
    Bounds3f FindAOBounds(float eps, Bounds3f modelBounds, Vec3f minArea, Vec3f maxArea);
    // Variant of FindAOBounds that lets you supply your own axial projected area bounds

    void CreateBitMaskFromTriangles
    (
        int             triCount,
        const int       indices[],
        const Vec3f     vertices[],
        const Bounds3f& bbox,
        int w, int h, int d,
        uint32_t        mask[]
    );
    ///< Initialises occupancy mask from the given mesh. 'indices' may be null.

    void CreateDirW8FromTriangles
    (
        int             triCount,
        const int       indices[],
        const Vec3f     vertices[],
        const Bounds3f& bbox,
        int w, int h, int d, 
        float           dirW8[]
    );
    ///< Initialises directional occlusion volumes from the given mesh. 'indices' may be null.

    struct cMesh
    {
        std::vector<Vec3f> mPositions;
        std::vector<Vec3f> mNormals;
        std::vector<Vec2f> mUVs;

        std::vector<int> mPositionIndices;
        std::vector<int> mNormalIndices;
        std::vector<int> mUVIndices;
    };

    bool ReadObjFile(FILE* file, cMesh* mesh);

    Bounds3f FindBounds(const cMesh& mesh);
}

#endif
