#include <vector>
#include <unordered_map>
#include <set>
#include "glm/glm.hpp"
#include "viewer.hpp"
#include <span>
#include <iostream>
#include <fstream>
#include <map>
#include <array>
#include <sstream>

#include <algorithm>

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
// using Vec3Span = std::span<Vec3>;
using IVec3List = std::vector<IVec3>;
using EdgeList = std::vector<Edge>;
// using EdgeSpan = std::span<Edge>;
using VecList = std::vector<Vec>;
using HalfEdgeList = std::vector<HalfEdge>;

// using IVec3Span = std::span<IVec3>;

class MeshHalfEdge
{

public:
    Vec3List vertexPos;    // n x 3 vertex positions
    Vec3List vertexNormal; // n x 3 vertex normals
    HalfEdgeList halfEdge; // e x 4 (PAIR, NEXT, HEAD, LEFT)

    IntList vertexHalfEdge; // Maps vertex index → one outgoing half-edge
    IntList FaceHalfEdge;   // Maps face index → one half-edge
    EdgeMap halfedgeMap;    // Edge lookup table

    // part1
    void addFace(const Face &face, IntList &prev);

    void handleBoundaryVertices(IntList &prev);

    void buildHalfEdgeStructure(VecList &faces);

    void triangulateFace(int faceIdx, IVec3List &triangleVertices, EdgeList &edges);

    void triangulateMesh(IVec3List &triangleVertices, EdgeList &edges);

    void sanity_check();

    void debugInfo(IVec3List &triangleVertices, EdgeList &edges);

    // part3
    void loadObjfile(const std::string &filename, Vec3List &vertices,

                     Vec3List &normals,
                     FaceList &faces);

    // part4
    Vec3 computeFaceNormal(int x, int y, int z);
    void computeVertexNormals();

    // part 5
    void smoothen_step(float λ, int maxVertex);
    void smoothen_laplacian(float λ, int iterations, int maxVertex = -1);
    void smoothen_taubin(float λ, float u, int iterations, int maxVertex = -1);

    void addNoise(const std::string &noiseType, float param);

    void computeEdgeStats(EdgeList &edges)
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

    // part 6
    // Extrudes a face outward by duplicating its vertices and adjusting connectivity.
    void extrudeFace(int faceidx, float multiplier = 1.0f);

    // part 8

    // (Only Geometrical) Moves a single vertex along its normal direction.
    void extrudeVertex(int v, float factor = 1.0f);

    // (Only Geometrical) Moves all neighboring vertices of `v` outward along `v`'s normal.
    void extrudeNeighbors(int v, float factor = 1.0f);

    // Duplicates and extrudes the neighbors of `v`, creating a new shell around `v`.
    void extrudeCopyNeighbors(int v, float factor = 1.0f);

    // (Only Geometrical) Moves `v` to the average position of its neighbors to create a local flattening effect.
    void flatten(int v);

    // (Only Geometrical) Adjusts neighboring vertices of `v` to lie on a plane perpendicular to its normal.
    void planarizeNeighbors(int v);
};

// // part2

// Grid along any axis at position X=a or Z=a etc
void generateCustomGrid(int m, int n, float a, int axis, const std::string &filename);

void generateGrid(int m, int n, const std::string &filename);
void generateSphere(int m, int n, const std::string &filename, int axis = 2, int direction = 1);
void generateCube(int m, int n, int o, const std::string &filename);
