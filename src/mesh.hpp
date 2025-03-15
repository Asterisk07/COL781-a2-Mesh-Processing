#include <vector>
#include <unordered_map>
#include <set>
#include "glm/glm.hpp"
#include "viewer.hpp"
#include <span>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>

#include <algorithm> // Required for std::reverse

inline void suppressStdout()
{
    std::cout.setstate(std::ios_base::failbit); // Disable cout
}

inline void restoreStdout()
{
    std::cout.clear(); // Restore cout
}

#include <cmath>

// Base case for recursion: print the last argument
template <typename T>
void print(const T &last)
{
    std::cout << last << std::endl;
}

template <typename T>
void printerr(const T &last)
{
    std::cerr << last << std::endl;
}

// Variadic template function to handle multiple arguments
template <typename T, typename... Args>
void print(const T &first, const Args &...rest)
{
    std::cout << first << " "; // Print the first argument with a space
    print(rest...);            // Recursively call print for remaining arguments
}

// Variadic template function to handle multiple arguments
template <typename T, typename... Args>
void printerr(const T &first, const Args &...rest)
{
    std::cerr << first << " "; // Print the first argument with a space
    printerr(rest...);         // Recursively call print for remaining arguments
}

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

using IntMap = std::unordered_map<int, int>;
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
    Vec3List vertexPos; // n x 3 vertex positions
    // Vec3List vertexPosFromFile;    // n x 3 vertex positions
    Vec3List vertexNormal; // n x 3 vertex normals
    // Vec3List vertexNormalFromFile; // n x 3 vertex normals
    // Vec3List faceNormals;  // m x 3 face normals
    HalfEdgeList halfEdge; // e x 4 (PAIR, NEXT, HEAD, LEFT)

    IntList vertexHalfEdge;     // Maps vertex index → one outgoing half-edge
    IntList FaceHalfEdge;       // Maps face index → one half-edge
    EdgeMap halfedgeMap;        // Edge lookup table
    IVec3List triangleVertices; // Explicit triangle indices for rendering
    EdgeList edges;             // Explicit triangle edges for rendering
    // EdgeSpan edges;             // Explicit triangle edges for rendering
    // VecList faceVertices; // Explicit face indices for rendering

    // part1
    void addFace(const Face &face, IntList &prev);
    void buildHalfEdgeStructure(IVec3Span triangles);
    void triangulateFace(int faceIdx);

    // void  resize(int n);
    void buildHalfEdgeStructure(VecList &faces);
    void initArray()
    {
        vertexHalfEdge.resize(vertexPos.size(), -1);
    }

    void sanity_check();
    void debugInfo()
    {

        std::cerr << " | Vertices " << vertexPos.size() << " | triangles " << triangleVertices.size() << " | Edges " << edges.size() << std::endl;

        // Print vertex index, position, and normal on the same line
        std::cerr << "Vertices" << std::endl;
        for (size_t idx = 0; idx < vertexPos.size(); idx++)
        {
            std::cerr << idx << ": (" << vertexPos[idx].x << ", " << vertexPos[idx].y << ", " << vertexPos[idx].z << ")"
                      << " | Normal: (" << vertexNormal[idx].x << ", " << vertexNormal[idx].y << ", " << vertexNormal[idx].z << ")"
                      << std::endl;
        }

        // Print face index, vertex indices, and face normal
        std::cerr << "Faces" << std::endl;

        int numFaces = FaceHalfEdge.size();
        for (int i = 0; i < numFaces; i++)
        {
            // triangulateFace(i);
            std::cerr << i << " :  " << FaceHalfEdge[i] << std::endl;
        }
        // FaceHalfEdge
        // for (size_t idx = 0; idx < triangleVertices.size(); idx++)
        // {
        //     std::cerr << idx << ": [" << triangleVertices[idx].x << ", " << triangleVertices[idx].y << ", " << triangleVertices[idx].z << "]"

        //               << std::endl;
        // }

        // count appeareances
        // {
        //     std::unordered_map<int, int> vertexCount;

        //     // Count vertex occurrences
        //     for (const auto &face : triangleVertices)
        //     {
        //         vertexCount[face.x]++;
        //         vertexCount[face.y]++;
        //         vertexCount[face.z]++;
        //     }

        //     // Print the results
        //     std::cerr << "Vertex appearances in triangles:\n";
        //     for (const auto &[vertex, count] : vertexCount)
        //     {
        //         std::cerr << "Vertex " << vertex << " appears in " << count << " triangles.\n";
        //     }
        // }

        // std::cout << "Vertices" << std::endl;
        // for (auto &i : vertexPos)
        // {
        //     std::cout << i.x << " " << i.y << " " << i.z << std::endl;
        // }
        // std::cout << "Normals" << std::endl;
        // for (auto &i : vertexNormal)
        // {
        //     std::cout << i.x << " " << i.y << " " << i.z << std::endl;
        // }

        // std::cout << "Face Normals" << std::endl;
        // for (auto &i : faceNormals)
        // {
        //     std::cout << i.x << " " << i.y << " " << i.z << std::endl;
        // }
        // std::cout << "Half Edges\n"
        //   << "PAIR NEXT HEAD LEFT" << std::endl;
        int flag = 1;

        int h;
        int f;
        // {
        //     int v = 7;
        //     h = vertexHalfEdge[v];
        //     // h = halfEdge[h][PAIR];
        //     // HalfEdge *h = v->halfEdge;
        //     printerr("Face traversal from vertex", v);
        //     do
        //     {

        //         // do something with h->left;
        //         f = halfEdge[h][LEFT];
        //         print("he ", h, "Face ", f);
        //         h = halfEdge[h][NEXT];
        //         // print("Next he ", h, "Face ", f, "head ", halfEdge[h][HEAD];);
        //         h = halfEdge[h][PAIR];

        //     } while (h != vertexHalfEdge[v]);
        // }
        // for (int i = 0; i < vertexHalfEdge.size(); i++)
        // {
        //     h = vertexHalfEdge[i];
        // std::cout << "Vertex " << i << " Points to " << h << "of head" << halfEdge[h][HEAD] << std::endl;
        // }

        // for (auto &i : halfEdge)
        printerr("Half edges : \nPAIR NEXT HEAD LEFT");
        for (int h = 0; h < halfEdge.size(); h++)
        {
            auto i = halfEdge[h];
            // if (h == 21 | h == 8 | h == 7)
            std::cerr << h << ":" << i[0] << " " << i[1] << " " << i[2] << " " << i[3] << std::endl;
            // if (i[0] == -1)
            // {
            //     std::cerr << "ERROR " << std::endl;
            //     flag = 0;
            //     break;
            // }
            // if (h==)
        }
        // if (flag)
        //     std::cout << "All half edge have pair" << std::endl;
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
    void handleBoundaryVertices(IntList &prev);
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
    void computeTriangleNormal(IVec3 &tri);
    Vec3 computeFaceNormal(int x, int y, int z);
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

    void extrudeFace(int faceidx, float multiplier = 1.0f);
};
// void findBoundaryEdges();
// // void extractEdgesFromFaces();

// // part2
void generateCustomGrid(int m, int n, float a, int axis, const std::string &filename);
void generateGrid(int m, int n, const std::string &filename);
void generateSphere(int m, int n, const std::string &filename);
void generateCube(int m, int n, int o, const std::string &filename);

// // part3
// void loadOBJ(std::string filename);
// void saveOBJ(std::string filename);

// auto resizeByN = [&](auto &vec)
// { vec.resize(vec.size() + n); };