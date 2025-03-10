#include "viewer.hpp"
#include "mesh.hpp"

namespace V = COL781::Viewer;
// using namespace glm;
using vec3 = Vec3;
using ivec3 = IVec3;
// using ivec2 = IVec2;
int main()
{

    vec3 vertices[] = {
        vec3(0.0, 0.5, 0.0),   // Top
        vec3(-0.5, 0.2, 0.0),  // Upper-left
        vec3(-0.3, -0.5, 0.0), // Bottom-left
        vec3(0.3, -0.5, 0.0),  // Bottom-right
        vec3(0.5, 0.2, 0.0)    // Upper-right
    };

    vec3 normals[] = {
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0),
        vec3(0.0, 0.0, 1.0)};

    // Define the face using the 5 vertices
    FaceList faces = {
        {0, 1, 2, 3, 4} // Pentagon face
    };

    // Create mesh using half-edge structure
    MeshHalfEdge mesh;

    mesh.vertexPos = vertices;

    mesh.vertexNormal = normals;

    // Build half-edge structure
    mesh.buildHalfEdgeStructure(faces);
    mesh.triangulateMesh();

    // ivec2 edges[] = {
    //     ivec2(0, 1),
    //     ivec2(1, 3),
    //     ivec2(3, 2),
    //     ivec2(2, 0),

    //     ivec2(4, 5),
    //     ivec2(5, 7),
    //     ivec2(7, 6),
    //     ivec2(6, 4),

    //     ivec2(0, 4),
    //     ivec2(1, 5),
    //     ivec2(2, 6),
    //     ivec2(3, 7)};
    // edges = 1;

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
