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
    suppressStdout();
    // const std::string filename = "meshes/cube.obj";
    // const std::string filename = "meshes/cube.obj";
    // const std::string filename = "meshes/try_tri.obj";
    const std::string filename = "meshes/try_new.obj";
    // const std::string filename = "meshes/try_square.obj";
    // const std::string filename = "meshes/try_pyramid.obj";
    // const std::string filename = "meshes/spot_control_mesh.obj";
    // const std::string filename = "meshes/bunny_1k.obj";
    // const std::string filename = "meshes/try_spot.obj";
    // generateGrid(3, 3, filename);
    // generateCustomGrid(4, 2, -1, 1, filename);
    int m = 7, n = 9;
    // int m = 6, n = 6;
    int axis = 1; // y axis
    int direction = -1;
    generateSphere(m, n, filename, axis, direction);

    // generateSphere(7, 9, filename);
    // generateCube(2, 3, 4, filename);
    // generateCube(3, 3, 3, filename);
    // generateCube(1, 3, 1, filename);
    // generateCube(4, 2, 5, filename);
    // generateCube(1, 1, 1, filename);
    mesh.loadObjfile(filename, vertices, texCoords, normals, faces);

    // // set Attribs
    // mesh.vertexPos = vertices;
    // mesh.vertexNormal = normals;

    // Build half-edge structure
    mesh.buildHalfEdgeStructure(faces);

    // Extrude cube
    // for (int i = 0; i < 4; i++)
    // {
    // mesh.extrudeFace(1, 1);
    // }

    mesh.sanity_check();
    std::cerr << "Here 0" << std::endl;
    // mesh.debugInfo();

    std::cerr << "Here" << std::endl;
    IVec3List triangleVertices; // Explicit triangle indices for rendering
    EdgeList edges;             // Explicit triangle edges for rendering
    mesh.triangulateMesh(triangleVertices, edges);
    std::cerr << "Here 2\n";
    // restoreStdout();
    // mesh.smoothen(0.8, 4);
    // mesh.computeEdgeStats(edges);
    // mesh.smoothen_taubin(0.4, -0.2, 50);
    // mesh.addNoise("gaussian", 0.05);
    // mesh.smoothen_laplacian(0.2, 8); // tk uncomment to smoothen
    mesh.smoothen_laplacian(-0.2, 15); // tk uncomment to smoothen
    // mesh.smoothen_laplacian(-0.1, 10);
    // mesh.extrudeVertex(0);
    mesh.extrudeVertex(-1, -0.5);
    mesh.extrudeVertex(0, 0.5);
    mesh.extrudeNeighbors(0, 0.5);
    // mesh.flatten(-1);
    mesh.planarizeNeighbors(-1);
    // mesh.extrudeVertex(-1, -0.5);
    // mesh.smoothen_laplacian(0.1, 10); // tk uncomment to smoothen

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
    v.setMesh(mesh.vertexPos.size(), triangleVertices.size(), edges.size(),
              mesh.vertexPos.data(), triangleVertices.data(), edges.data(), mesh.vertexNormal.data());

    v.view();
}
