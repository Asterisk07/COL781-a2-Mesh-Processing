#include "mesh.hpp" // ✅ Always include the corresponding header

// A unit square in the xy-plane, divided into a grid of m rows and n columns. Each grid cell will be a 1/m × 1/n rectangle, for a total of (m+1)(n+1) vertices and mn quad faces.
// void generateGridMesh(int m, int n) {

// };

// void loadObjfile(){

// }

// Base case for recursion: print the last argument
template <typename T>
void print(const T &last)
{
    std::cout << last << std::endl;
}

// Variadic template function to handle multiple arguments
template <typename T, typename... Args>
void print(const T &first, const Args &...rest)
{
    std::cout << first << " "; // Print the first argument with a space
    print(rest...);            // Recursively call print for remaining arguments
}

// No need to include <vector>, <unordered_map>, etc. again if they are in mesh.hpp
// part1
Key enkey(int x, int y)
{
    return ((Key(x) << 32) | y);
}

// void MeshHalfEdge::resize(int m, int n)
// {

//     // int m = vertexPos.size();
//     vertexPos.resize(this->size() + n);
//     vertexNormal.resize(this->size() + n);

//     halfEdge.resize(this->size() + m);
// }

void MeshHalfEdge::addFace(const Face &face)
{
    int n = face.size();
    assert(n > 0);
    int key;
    // HalfEdge *h;

    int faceIdx = FaceHalfEdge.size();

    int m = halfEdge.size(); // Get starting index for new half-edges

    halfEdge.resize(m + n, {-1, -1, -1, -1}); // Resize for new half-edges
    // vertexHalfEdge.resize(vertexHalfEdge.size()+)
    int vertex;
    int nextVertex;
    int pair;
    int h;

    FaceHalfEdge.push_back(m);
    for (int i = 0; i < n; i++)
    {
        h = m + i;
        vertex = face[i];
        nextVertex = face[(i + 1) % n];
        std::cout << "vertex is " << vertex << std::endl;
        assert(vertex < vertexHalfEdge.size());
        assert(nextVertex < vertexHalfEdge.size());
        vertexHalfEdge[vertex] = h;

        key = enkey(nextVertex, vertex);
        if (halfedgeMap.find(key) != halfedgeMap.end())
        {
            pair = halfedgeMap[key];
            halfEdge[h][PAIR] = pair;
            halfEdge[pair][PAIR] = h;
        }
        else
        {
            key = enkey(vertex, nextVertex);
            halfedgeMap[key] = h;
        }

        // halfEdge[h][PREV] = (i - 1 + n) % n + m;

        halfEdge[h][NEXT] = (i + 1) % n + m;

        halfEdge[h][LEFT] = faceIdx;
        halfEdge[h][HEAD] = nextVertex;
    }
};

void MeshHalfEdge::buildHalfEdgeStructure(VecList &faces)
{
    initArray();
    for (Face &face : faces)
    {
        addFace(face);
    }
}

Face triangle_to_face(IVec3 &tri)
{
    return {tri.x, tri.y, tri.z}; // Directly initialize the vector
}

void MeshHalfEdge::buildHalfEdgeStructure(IVec3Span triangles)
{

    initArray();
    // convert to faces
    // Face face;
    for (auto &tri : triangles)
    {
        addFace(triangle_to_face(tri));
    }
    // buildHalfEdgeStructure(faces);
}

// void MeshHalfEdge::findBoundaryEdges();

void print()
{
}

void MeshHalfEdge::triangulateFace(int faceIdx)
{
    // triangulates face and adds triangle to triangleVertices
    // int h = FaceHalfEdge[faceIdx];
    // HalfEdge h =
    HalfEdge h;

    h = halfEdge[FaceHalfEdge[faceIdx]];
    int x = -1, y = -1, z = -1;
    x = h[HEAD];
    print("currently h[HEAD] = ", h[HEAD]);

    h = halfEdge[h[NEXT]];
    print("x,y,z", x, y, z);

    while (h[HEAD] != x)
    {
        z = h[HEAD];
        if (y == -1)
        {
            edges.push_back(Edge(x, z));
        }
        else
        {
            // edges.push_back({y, z});
            edges.push_back(Edge(y, z));

            // IVec3 tri(x, y, z);
            triangleVertices.push_back(IVec3(x, y, z));
        }

        y = z;
        h = halfEdge[h[NEXT]];
        print("x,y,z", x, y, z);
    }
}
// void MeshHalfEdge::extractEdgesFromFaces();

// // part2
// void MeshHalfEdge::generateGridMesh(int m, int n);
// void MeshHalfEdge::generateSphereMesh(int m, int n);
// void MeshHalfEdge::enerateCubeMesh(int m, int n, int o);

// // part3
// void MeshHalfEdge::loadOBJ(std::string filename);
// void MeshHalfEdge::saveOBJ(std::string filename);

// // part4
// void MeshHalfEdge::computeFaceNormals();
// void MeshHalfEdge::computeVertexNormals();
