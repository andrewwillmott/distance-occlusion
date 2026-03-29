//------------------------------------------------------------------------------
// Purpose: Quick and dirty mesh support
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#define  _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include "MeshSupport.hpp"

#include <assert.h>
#include <float.h>
#include <stdio.h>

#ifdef _MSC_VER
    #pragma warning (disable: 4244)
    #define strcasecmp _stricmp
    #define strtok_r strtok_s
#endif

using namespace MSL;

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

Bounds3 MSL::FindBounds(const cMesh& mesh)
{
    Bounds3 bbox = { { +FLT_MAX, +FLT_MAX, +FLT_MAX }, { -FLT_MAX, -FLT_MAX, -FLT_MAX } };

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


