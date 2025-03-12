#include "viewer.hpp"
#include "mesh.hpp"

namespace V = COL781::Viewer;
// using namespace glm;
using vec3 = Vec3;
using ivec3 = IVec3;
// using ivec2 = IVec2;
int main()
{

    // Define a single triangle
    // IVec3 faces[] = {
    //     {0, 1, 2}};

    // // Define vertex positions
    // Vec3 vertices[] = {
    //     Vec3(-0.5f, -0.5f, 0.0f), // v0
    //     Vec3(0.5f, -0.5f, 0.0f),  // v1
    //     Vec3(0.0f, 0.5f, 0.0f)    // v2
    // };

    // Vec3 normals[] = {
    //     Vec3(0.0, 0.0, -1.0),
    //     Vec3(0.0, 0.0, -1.0),
    //     Vec3(0.0, 0.0, -1.0),
    // };

    /* ------------------------------------ ------------------------------------*/

    Vec3List vertices;
    Vec2List texCoords;
    Vec3List normals;
    FaceList faces;

    // Create mesh using half-edge structure
    MeshHalfEdge mesh;
    // const std::string filename = "meshes/cube.obj";
    // const std::string filename = "meshes/cube.obj";
    // const std::string filename = "meshes/try_tri.obj";
    const std::string filename = "meshes/try_pyramid_regular.obj";
    // const std::string filename = "meshes/try_pyramid.obj";
    // const std::string filename = "meshes/spot_control_mesh.obj";
    // const std::string filename = "meshes/bunny_1k.obj";
    // const std::string filename = "meshes/try_spot.obj";
    mesh.loadObjfile(filename, vertices, texCoords, normals, faces);

    // // set Attribs
    // mesh.vertexPos = vertices;
    // mesh.vertexNormal = normals;

    // Build half-edge structure
    mesh.buildHalfEdgeStructure(faces);
    mesh.triangulateMesh();
    // mesh.smoothen(0.8, 4);
    mesh.computeEdgeStats();
    mesh.smoothen_taubin(0.2, -0.1, 5);
    // mesh.smoothen_laplacian(0.2, 5);
    mesh.computeEdgeStats();
    // mesh.computeVertexNormals();

    /* ------------------------------------ ------------------------------------*/

    // Initialize the viewer
    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480))
    {
        return EXIT_FAILURE;
    }

    // Set mesh data
    // mesh.debugInfo();
    v.setMesh(mesh.vertexPos.size(), mesh.triangleVertices.size(), mesh.edges.size(),
              mesh.vertexPos.data(), mesh.triangleVertices.data(), mesh.edges.data(), mesh.vertexNormal.data());

    v.view();
}
