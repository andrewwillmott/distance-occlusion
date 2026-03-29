//------------------------------------------------------------------------------
// Purpose: Quick and dirty mesh support
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#ifndef MESH_SUPPORT_H
#define MESH_SUPPORT_H

#include "VLMini.hpp"
#include <vector>
#include <stdio.h>

namespace MSL
{
    using DOL::Bounds3;

    // Mesh
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

    Bounds3 FindBounds(const cMesh& mesh);
}

#endif
