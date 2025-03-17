#include "viewer.hpp"
#include "mesh.hpp"

namespace V = COL781::Viewer;
// using namespace glm;
using vec3 = Vec3;
using ivec3 = IVec3;
// using ivec2 = IVec2;
int main()
{

    Vec3List vertices;
    Vec3List normals;
    FaceList faces;

    // Create mesh using half-edge structure
    MeshHalfEdge mesh;
    suppressStdout();

    const std::string filename = "meshes/try_new.obj";
    // const std::string filename = "meshes/bunny_1k.obj";
    // const std::string filename = "meshes/cube.obj";
    // const std::string filename = "meshes/spot_control_mesh.obj";

    // generateSphere(7, 9, filename);
    generateGrid(6, 4, filename);
    // generateCube(2, 3, 4, filename);

    // generateSphere(4, 2, filename);
    // generateGrid(1, 1, filename);
    // generateCube(3, 3, 3, filename);
    // generateCube(3, 3, 3, filename);
    // generateCube(4, 2, 5, filename);
    // generateCube(1, 1, 1, filename);
    mesh.loadObjfile(filename, vertices, normals, faces);

    // Build half-edge structure
    mesh.buildHalfEdgeStructure(faces);

    mesh.sanity_check();
    std::cerr << "Here 0" << std::endl;
    // mesh.debugInfo();

    std::cerr << "Here" << std::endl;
    IVec3List triangleVertices; // Explicit triangle indices for rendering
    EdgeList edges;             // Explicit triangle edges for rendering
    mesh.triangulateMesh(triangleVertices, edges);
    std::cerr << "Here 2\n";
    mesh.computeVertexNormals();

    /* ------------------------------------ ------------------------------------*/

    // Initialize the viewer
    V::Viewer v;
    if (!v.initialize("Mesh viewer", 640, 480))
    {
        return EXIT_FAILURE;
    }

    // Set mesh data
    mesh.debugInfo(triangleVertices, edges);
    v.setMesh(mesh.vertexPos.size(), triangleVertices.size(), edges.size(),
              mesh.vertexPos.data(), triangleVertices.data(), edges.data(), mesh.vertexNormal.data());

    v.view();
}
