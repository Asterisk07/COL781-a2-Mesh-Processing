#include <vector>
#include <unordered_map>
#include <set>
#include "glm/glm.hpp"
#include "viewer.hpp"
#include <span>
#include <iostream>
#include <fstream>
#include <sstream>

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
using FaceList = std::vector<Face>;

using Vec3 = glm::vec3;
using Vec2 = glm::vec2;
using IVec3 = glm::ivec3;
using Edge = glm::ivec2;
using Vec = std::vector<int>;
using HalfEdge = std::array<int, 4>;

using Vec3List = std::vector<Vec3>;
using Vec2List = std::vector<Vec2>;
using Vec3Span = std::span<Vec3>;
using IVec3List = std::vector<IVec3>;
// using IVec3Array = IVec3[];
using EdgeList = std::vector<Edge>;
using EdgeSpan = std::span<Edge>;
using VecList = std::vector<Vec>;
using HalfEdgeList = std::vector<HalfEdge>;

using IVec3Span = std::span<IVec3>;

class MeshHalfEdge
{

public:
    Vec3Span vertexPos; // n x 3 vertex positions
    // Vec3List vertexPosFromFile;    // n x 3 vertex positions
    Vec3Span vertexNormal; // n x 3 vertex normals
    // Vec3List vertexNormalFromFile; // n x 3 vertex normals
    Vec3List faceNormals;  // m x 3 face normals
    HalfEdgeList halfEdge; // e x 4 (PAIR, NEXT, HEAD, LEFT)

    // Optional but useful
    IntList vertexHalfEdge;     // Maps vertex index → one outgoing half-edge
    IntList FaceHalfEdge;       // Maps face index → one half-edge
    EdgeMap halfedgeMap;        // Edge lookup table
    IVec3List triangleVertices; // Explicit triangle indices for rendering
    EdgeList edges;             // Explicit triangle edges for rendering
    // EdgeSpan edges;             // Explicit triangle edges for rendering
    VecList faceVertices; // Explicit face indices for rendering

    // part1
    void addFace(const Face &face);
    void buildHalfEdgeStructure(IVec3Span triangles);
    void triangulateFace(int faceIdx);

    // void  resize(int n);
    void buildHalfEdgeStructure(VecList &faces);
    void initArray()
    {
        vertexHalfEdge.resize(vertexPos.size(), -1);
    }

    void debugInfo()
    {

        std::cout << " | Vertices " << vertexPos.size() << " | triangles " << triangleVertices.size() << " | Edges " << edges.size() << std::endl;

        std::cout << "Vertices" << std::endl;
        for (auto &i : vertexPos)
        {
            std::cout << i.x << " " << i.y << " " << i.z << std::endl;
        }
        std::cout << "Normals" << std::endl;
        for (auto &i : vertexNormal)
        {
            std::cout << i.x << " " << i.y << " " << i.z << std::endl;
        }
        std::cout << "Edges" << std::endl;
        for (auto &i : edges)
        {
            std::cout << i.x << " " << i.y << std::endl;
        }
        std::cout << "Triangles" << std::endl;
        for (auto &i : triangleVertices)
        {
            std::cout << i.x << " " << i.y << " " << i.z << std::endl;
        }
    }

    void triangulateMesh()
    {
        int numFaces = FaceHalfEdge.size();
        for (int i = 0; i < numFaces; i++)
        {
            triangulateFace(i);
        }
    }

    // void loadObjfile(const std::string &filename);

    void loadObjfile(const std::string &filename, Vec3List &vertices,
                     Vec2List &texCoords,
                     Vec3List &normals,
                     FaceList &faces);
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
