#include <vector>
#include <unordered_map>
#include <set>
#include "glm/glm.hpp"
#include "viewer.hpp"
#include <span>
#include <iostream>
#include <fstream>
#include <sstream>
// #include <glm/gtc/magnitude.hpp> // Optional, part of GLM

inline double edgeLength(const glm::vec3 &v1, const glm::vec3 &v2)
{
    return glm::length(v1 - v2); // Corrected for GLM
}

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

        std::cout << "Face Normals" << std::endl;
        for (auto &i : faceNormals)
        {
            std::cout << i.x << " " << i.y << " " << i.z << std::endl;
        }
        std::cout << "Half Edges\n"
                  << "PAIR NEXT HEAD LEFT" << std::endl;
        for (auto &i : halfEdge)
        {
            std::cout << i[0] << " " << i[1] << " " << i[2] << " " << i[3] << std::endl;
            if (i[0] == -1)
            {
                std::cout << "ERROR " << std::endl;
                break;
            }
        }
        // std::cout << "Edges" << std::endl;
        // for (auto &i : edges)
        // {
        //     std::cout << i.x << " " << i.y << std::endl;
        // }
        // std::cout << "Triangles" << std::endl;
        // for (auto &i : triangleVertices)
        // {
        //     std::cout << i.x << " " << i.y << " " << i.z << std::endl;
        // }
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

    // part3
    void loadObjfile(const std::string &filename, Vec3List &vertices,
                     Vec2List &texCoords,
                     Vec3List &normals,
                     FaceList &faces);

    // part4
    void computeFaceNormal(IVec3 &tri);
    void computeVertexNormals(); // tk incomplete : handle boundaries

    // part 5
    void smoothen_laplacian(float λ, int iterations);
    void smoothen_taubin(float λ, float u, int iterations);
    void smoothen_step(float λ);

    void computeEdgeStats()
    {
        std::vector<double> lengths;
        double sum = 0.0, sumSq = 0.0;

        for (auto &e : edges)
        {
            double len = edgeLength(vertexPos[e.x], vertexPos[e.y]);
            lengths.push_back(len);
            sum += len;
            sumSq += len * len;
        }

        double mean = sum / lengths.size();
        double variance = (sumSq / lengths.size()) - (mean * mean);
        double stddev = sqrt(variance);

        std::cout << "Mean Edge Length: " << mean << std::endl;
        std::cout << "Standard Deviation: " << stddev << std::endl;
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

// auto resizeByN = [&](auto &vec)
// { vec.resize(vec.size() + n); };