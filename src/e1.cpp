#include "viewer.hpp"
#include "mesh.hpp"

namespace V = COL781::Viewer;
// using namespace glm;
using vec3 = Vec3;
using ivec3 = IVec3;
// using ivec2 = IVec2;
int main()
{
    // Create mesh using half-edge structure
    MeshHalfEdge mesh;

    // Define a single triangle
    IVec3List triangles = {
        {0, 1, 2}};

    // Define vertex positions
    Vec3 vertices[] = {
        Vec3(-0.5f, -0.5f, 0.0f), // v0
        Vec3(0.5f, -0.5f, 0.0f),  // v1
        Vec3(0.0f, 0.5f, 0.0f)    // v2
    };

    Vec3 normals[] = {
        Vec3(0.0, 0.0, -1.0),
        Vec3(0.0, 0.0, -1.0),
        Vec3(0.0, 0.0, -1.0),
    };

    mesh.vertexPos = vertices;

    mesh.vertexNormal = normals;

    // Build half-edge structure
    mesh.buildHalfEdgeStructure(triangles);
    mesh.triangulateMesh();

    // Extract mesh data for rendering
    // Vec3L vertices = mesh.triangleVertices;

    // IVec3List triList = mesh.triangleVertices;

    // Extract edges for wireframe rendering
    // mesh.extractEdgesFromFaces();
    // std::vector<ivec2> edges = mesh.edges;

    // Compute normals
    // mesh.computeFaceNormals();
    // mesh.computeVertexNormals();
    // std::vector<vec3> normals = mesh.vertexNormal;

    // Initialize the viewer
    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480))
    {
        return EXIT_FAILURE;
    }

    // Set mesh data
    mesh.debugInfo();
    v.setMesh(mesh.vertexPos.size(), mesh.triangleVertices.size(), mesh.edges.size(),
              mesh.vertexPos.data(), mesh.triangleVertices.data(), mesh.edges.data(), mesh.vertexNormal.data());

    v.view();
}
