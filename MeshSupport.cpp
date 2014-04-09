//------------------------------------------------------------------------------
// Purpose: Quick and dirty mesh support
// Author:  Andrew Willmott
//------------------------------------------------------------------------------

#include "MeshSupport.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

using namespace MSL;

namespace
{
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

    inline float len(Vec3f v)
    {
        return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    Vec2f SquareToTriangle(Vec2f c)
    {
        float t = sqrt(c.y);

        Vec2f result =
        {
            t * c.x,
            1.0f - t
        };

        return result;
    }

    struct cHalton2 : public Vec2f
    {
        uint32_t mBase2;
        uint32_t mBase3;
        
        cHalton2() :
            mBase2(0),
            mBase3(0)
        {
            x = 0.0f;
            y = 0.0f;
        }
        
        int operator++()
        {
            uint32_t mask = 0x1;
            float s = 0.5f;
            mBase2++;

            while ((mBase2 & mask) == 0)
            {
                mask = mask << 1;

                x -= s;
                s *= 0.5f;
            }

            x += s;

            mask = 0x3;
            uint32_t add  = 0x1;
            s = 1.0f / 3.0f;
            mBase3++;

            while ((mBase3 & mask) == mask)
            {
                mBase3 += add;          // force carry into next 2-bit digit

                mask = mask << 2;
                add  = add  << 2;
                
                y -= 2.0f * s;
                s *= 1.0f / 3.0f;
            }

            y += s;

            return mBase2; // return the index of this sequence point
        }
    };

    const int kDx[8] = { +1, -1, +1, -1, +1, -1, +1, -1 };
    const int kDy[8] = { +1, +1, -1, -1, +1, +1, -1, -1 };
    const int kDz[8] = { +1, +1, +1, +1, -1, -1, -1, -1 };

    int RasteriseTriangleToDirW
    (
        const Vec3f&    v0,
        const Vec3f&    v1,
        const Vec3f&    v2,
        int w, int h, int d, 
        float           dirW[]
    )
    {
        int rowStride = w;
        int sliceStride = rowStride * h;
        int volumeStride = sliceStride * d;
        int dirWSize = volumeStride * 8;

        Vec3f an = TriDoubleAreaNormal(v0, v1, v2);

        float area = 0.5f * len(an);

        int samples = CeilToInt32(area) * 4;

        int dirOffset[8];
        float dirValue[8];
        int dirCount = 0;

        // We have two conflicting desires here -- we want the reconstruction of the occlusion direction
        // to be the same as the triangle normal, and ideally the occlusion amount = 0.5 = the hemisphere.
        // Luckily for our basis M, MtM = 8I, and more so, if you zero negative weights,
        // MtM = 4I. As we treat each weight as an octant coverage value however, we have an issue
        // in that e_z -> b_i = 1, and e_111 -> b_i = (3, 1, 1, 1) / sqrt(3), whose max coverage > 1, which
        // our sweep algorithm doesn't handle.

        float bs = 0.577350269f / area; //  1 / ||a|| sqrt(3)

        for (int i = 0; i < 8; i++)
        {
            float bc = kDx[i] * an.x + kDy[i] * an.y + kDz[i] * an.z;
            
            if (bc > 1e-6f)
            {
                dirOffset[dirCount] = i * volumeStride;
                dirValue [dirCount] = bc * bs;
                dirCount++;
            }
        }

        // This is a total goof. Rather than write a chunk of rasterisation code, just use the halton sequence
        // to sample the triangle and splat those points into the volume.
        cHalton2 hs;

        for (int i = 0; i < samples; i++)
        {
            ++hs;
            Vec2f ths = SquareToTriangle(hs);
            float s = ths.x;
            float t = ths.y;
            float u = 1 - s - t;

            Vec3f p =
            {
                s * v0.x + t * v1.x + u * v2.x,
                s * v0.y + t * v1.y + u * v2.y,
                s * v0.z + t * v1.z + u * v2.z
            };

            int x = FloorToInt32(p.x);
            int y = FloorToInt32(p.y);
            int z = FloorToInt32(p.z);

            if (x < 0 || x >= w)
                continue;
            if (y < 0 || y >= h)
                continue;
            if (z < 0 || z >= d)
                continue;

            int index = z * sliceStride + y * rowStride + x;
            
            for (int id = 0; id < dirCount; id++)
            {
                int dirIndex = dirOffset[id] + index;

                assert(dirIndex >= 0 && dirIndex < dirWSize);

                if (dirW[dirIndex] < dirValue[id])
                    dirW[dirIndex] = dirValue[id];
            }
        }

        return samples;
    }
}

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


bool MSL::CreateDirW8FromTriangles
(
    int             triCount,
    const Vec3f     vertices[],
    const Bounds3f& bbox,
    int w, int h, int d, 
    float           dirW8[]
)
{
    Vec3f invWidth =
    {
        float(w) / (bbox.mMax.x - bbox.mMin.x),
        float(h) / (bbox.mMax.y - bbox.mMin.y),
        float(d) / (bbox.mMax.z - bbox.mMin.z)
    };

    for (int i = 0, iv = 0; i < triCount; i++, iv += 3)
    {
        Vec3f v0 = vertices[iv + 0];
        Vec3f v1 = vertices[iv + 1];
        Vec3f v2 = vertices[iv + 2];

        v0.x = (v0.x - bbox.mMin.x) * invWidth.x;
        v0.y = (v0.y - bbox.mMin.y) * invWidth.y;
        v0.z = (v0.z - bbox.mMin.z) * invWidth.z;

        v1.x = (v1.x - bbox.mMin.x) * invWidth.x;
        v1.y = (v1.y - bbox.mMin.y) * invWidth.y;
        v1.z = (v1.z - bbox.mMin.z) * invWidth.z;

        v2.x = (v2.x - bbox.mMin.x) * invWidth.x;
        v2.y = (v2.y - bbox.mMin.y) * invWidth.y;
        v2.z = (v2.z - bbox.mMin.z) * invWidth.z;

        RasteriseTriangleToDirW
        (
            v0, v1, v2,
            w, h, d, dirW8
        );
    }

    return true;
}


bool MSL::CreateDirW8FromTriangles
(
    int             triCount,
    const int       vertexIndices[],
    const Vec3f     vertices[],
    const Bounds3f& bbox,
    int w, int h, int d, 
    float           dirW8[]
)
{
    Vec3f invWidth =
    {
        float(w) / (bbox.mMax.x - bbox.mMin.x),
        float(h) / (bbox.mMax.y - bbox.mMin.y),
        float(d) / (bbox.mMax.z - bbox.mMin.z)
    };

    for (int i = 0, iv = 0; i < triCount; i++, iv += 3)
    {
        int vi0 = vertexIndices[iv + 0];
        int vi1 = vertexIndices[iv + 1];
        int vi2 = vertexIndices[iv + 2];

        Vec3f v0 = vertices[vi0];
        Vec3f v1 = vertices[vi1];
        Vec3f v2 = vertices[vi2];

        v0.x = (v0.x - bbox.mMin.x) * invWidth.x;
        v0.y = (v0.y - bbox.mMin.y) * invWidth.y;
        v0.z = (v0.z - bbox.mMin.z) * invWidth.z;

        v1.x = (v1.x - bbox.mMin.x) * invWidth.x;
        v1.y = (v1.y - bbox.mMin.y) * invWidth.y;
        v1.z = (v1.z - bbox.mMin.z) * invWidth.z;

        v2.x = (v2.x - bbox.mMin.x) * invWidth.x;
        v2.y = (v2.y - bbox.mMin.y) * invWidth.y;
        v2.z = (v2.z - bbox.mMin.z) * invWidth.z;

        RasteriseTriangleToDirW
        (
            v0, v1, v2,
            w, h, d,
            dirW8
        );
    }

    return true;
}

namespace
{
    void InPlaceTriangulate(int numVerts, std::vector<int>& indices)
    {
        // Assume polygon of numVerts indices at the end of 'indices', and triangulate it in place.

        int baseVertexIndex = indices.size() - numVerts;
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

    bool ObjectCommand(cMesh* mesh, int argc, const char* va[])
    {
        return true;
    }

    bool GroupCommand(cMesh* mesh, int argc, const char* va[])
    {
        return true;
    }

    bool SmoothingGroupCommand(cMesh* mesh, int argc, const char* va[])
    {
        return true;
    }

    bool MaterialCommand(cMesh* mesh, int argc, const char* va[])
    {
        return true;
    }
    bool MaterialLibraryCommand(cMesh* mesh, int argc, const char* va[])
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
        maxArgs--;  // always reserve the last spot for 0 terminator
        
        for (int argc = 0; argc < maxArgs; argc++)
        {
            argv[argc] = strsep(&buffer, " \t\n\r");

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

    char* lineBuffer = 0;
    size_t lineBufferSize = 0;
    ssize_t lineSize;

    while ((lineSize = getline(&lineBuffer, &lineBufferSize, file)) > 0)
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
