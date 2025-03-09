#include <vector>
#include <unordered_map>
#include <set>
#include "glm/glm.hpp"
#include "viewer.hpp"
#include <span>
#include <iostream>

enum HalfEdgeIndices
{
    PAIR = 0,
    NEXT,
    HEAD,
    LEFT
};

using Key = uint64_t;
using EdgeMap = std::unordered_map<Key, int>;

using Face = std::vector<int>;
using IntList = std::vector<int>;

using Vec3 = glm::vec3;
using IVec3 = glm::ivec3;
using Edge = glm::ivec2;
using Vec = std::vector<int>;
using HalfEdge = std::array<int, 4>;

using Vec3List = std::vector<Vec3>;
using Vec3Span = std::span<Vec3>;
using IVec3List = std::vector<IVec3>;
using EdgeList = std::vector<Edge>;
using VecList = std::vector<Vec>;
using HalfEdgeList = std::vector<HalfEdge>;

class MeshHalfEdge
{

public:
    Vec3Span vertexPos;    // n x 3 vertex positions
    Vec3Span vertexNormal; // n x 3 vertex normals
    // Vec3List vertexNormal; // n x 3 vertex normals
    Vec3List faceNormals;  // m x 3 face normals
    HalfEdgeList halfEdge; // e x 4 (PAIR, NEXT, HEAD, LEFT)

    // Optional but useful
    IntList vertexHalfEdge;     // Maps vertex index → one outgoing half-edge
    IntList FaceHalfEdge;       // Maps face index → one half-edge
    EdgeMap halfedgeMap;        // Edge lookup table
    IVec3List triangleVertices; // Explicit triangle indices for rendering
    EdgeList edges;             // Explicit triangle edges for rendering
    VecList faceVertices;       // Explicit face indices for rendering

    // part1
    void addFace(const Face &face);
    void buildHalfEdgeStructure(IVec3List &triangles);
    void triangulateFace(int faceIdx);

    // void  resize(int n);
    void buildHalfEdgeStructure(VecList &faces);
    void initArray()
    {
        vertexHalfEdge.resize(vertexPos.size(), -1);
    }

    void debugInfo()
    {
        std::cout << vertexPos.size() << triangleVertices.size() << edges.size() << std::endl;
    }

    void triangulateMesh()
    {
        int numFaces = FaceHalfEdge.size();
        for (int i = 0; i < numFaces; i++)
        {
            triangulateFace(i);
        }
    }
};
// void findBoundaryEdges();
// // void extractEdgesFromFaces();

// // part2
// void generateGridMesh(int m, int n);
// void generateSphereMesh(int m, int n);
// void enerateCubeMesh(int m, int n, int o);

// // part3
// void loadOBJ(std::string filename);
// void saveOBJ(std::string filename);

// // part4
// void computeFaceNormals();
// void computeVertexNormals();

// auto resizeByN = [&](auto &vec)
// { vec.resize(vec.size() + n); };
