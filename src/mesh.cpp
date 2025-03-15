#include "mesh.hpp" // ✅ Always include the corresponding header

// A unit square in the xy-plane, divided into a grid of m rows and n columns. Each grid cell will be a 1/m × 1/n rectangle, for a total of (m+1)(n+1) vertices and mn quad faces.
// void generateGridMesh(int m, int n) {

// };

void MeshHalfEdge::loadObjfile(const std::string &filename, Vec3List &vertices,
                               Vec2List &texCoords,
                               Vec3List &normals,
                               FaceList &faces)
{
    // Vec3List vertices;
    // Vec2List texCoords;
    // Vec3List normals;
    // FaceList faces;

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Cannot open OBJ file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "v")
        { // Vertex positions
            Vec3 v;
            iss >> v.x >> v.y >> v.z;
            std::cout << "\tOne row of  vertices : " << v.x << " " << v.y << " " << v.z << std::endl;
            vertices.push_back(v);

            // std::cerr << "Current Vertices" << std::endl;
            // for (auto &i : vertices)
            // {
            // std::cerr << i.x << " " << i.y << " " << i.z << std::endl;
            // }
        }
        else if (prefix == "vt")
        { // Texture coordinates
            Vec2 vt;
            iss >> vt.x >> vt.y;
            texCoords.push_back(vt);
        }
        else if (prefix == "vn")
        { // Vertex normals
            Vec3 vn;
            iss >> vn.x >> vn.y >> vn.z;
            normals.push_back(vn);
        }
        else if (prefix == "f")
        { // Faces
            Face face;
            std::string vertexData;
            while (iss >> vertexData)
            {
                std::istringstream vss(vertexData);
                std::string vIdx, vtIdx, vnIdx;
                std::getline(vss, vIdx, '/');  // Read vertex index
                std::getline(vss, vtIdx, '/'); // Read texture index (if exists)
                std::getline(vss, vnIdx, '/'); // Read normal index (if exists)

                // face.vertexIndices.push_back(std::stoi(vIdx) - 1); // OBJ indices start at 1
                face.push_back(std::stoi(vIdx) - 1); // OBJ indices start at 1

                // if (!vtIdx.empty()) face.texIndices.push_back(std::stoi(vtIdx) - 1);
                // if (!vnIdx.empty())
                //     face.normalIndices.push_back(std::stoi(vnIdx) - 1);
            }
            faces.push_back(face);
        }
    }

    file.close();
    this->vertexPos = std::move(vertices);
    // this->vertexPos = vertices;
    // // mesh.texCoords = texCoords;
    if (normals.size() < vertices.size())
    {
        normals.resize(vertices.size(), Vec3(0.0f, 0.0f, 0.0f)); // Extend with zero normals
    }
    this->vertexNormal = std::move(normals);
    // this->vertexNormal = normals;

    // this->buildHalfEdgeStructure(faces);
    // this->triangulateMesh();
}

// // Base case for recursion: print the last argument
// template <typename T>
// void print(const T &last)
// {
//     std::cout << last << std::endl;
// }

// // Variadic template function to handle multiple arguments
// template <typename T, typename... Args>
// void print(const T &first, const Args &...rest)
// {
//     std::cout << first << " "; // Print the first argument with a space
//     print(rest...);            // Recursively call print for remaining arguments
// }

// No need to include <vector>, <unordered_map>, etc. again if they are in mesh.hpp
// part1
Key enkey(int x, int y)
{
    // print("KEY : ", x, " + ", y, " = ", (Key(x) << 32) | y);

    // return ((Key(x) * 10) + y);
    return ((Key(x) << 32) | y);
}

// void MeshHalfEdge::resize(int m, int n)
// {

//     // int m = vertexPos.size();
//     vertexPos.resize(this->size() + n);
//     vertexNormal.resize(this->size() + n);

//     halfEdge.resize(this->size() + m);
// }

// void MeshHalfEdge::extrudeFace(int faceidx)
// {
// }
void MeshHalfEdge::extrudeFace(int faceidx, float multiplier)
{

    int h = FaceHalfEdge[faceidx];

    int u = -1, v = -1, w = -1;
    u = halfEdge[h][HEAD];
    h = halfEdge[h][NEXT];
    v = halfEdge[h][HEAD];
    h = halfEdge[h][NEXT];

    int H_old = halfEdge.size();
    FaceHalfEdge[faceidx] = H_old + 3;
    int V_old = vertexPos.size();

    int start = h;
    int x;

    do
    {

        w = halfEdge[h][HEAD];

        int V = vertexPos.size();
        int H = halfEdge.size();
        int F = FaceHalfEdge.size();
        FaceHalfEdge.push_back(H);

        // halfEdge[h][NEXT] = H;

        vertexPos.push_back(vertexPos[v]);
        vertexNormal.push_back(vertexNormal[v]);

        printerr("Now size of vertices : ", vertexPos.size());

        vertexHalfEdge.push_back(H + 3);
        halfEdge.push_back(std::array<int, 4>{H + 6, H + 1, V + 1, F});
        halfEdge.push_back(std::array<int, 4>{H + 3, H + 2, V, F});
        halfEdge.push_back(std::array<int, 4>{H - 4, h, v, F});
        halfEdge.push_back(std::array<int, 4>{H + 1, H + 7, V + 1, faceidx});

        if (halfEdge[h][NEXT] == start)
        {
            // i.e. this is last one
            halfEdge[H][PAIR] = H_old + 2;
            halfEdge[H_old + 2][PAIR] = H;

            halfEdge[H][HEAD] = V_old;
            halfEdge[H + 3][HEAD] = V_old;

            halfEdge[H + 3][NEXT] = H_old + 3;
        }

        Vec3 normal = {0.0f, 0.0f, 0.0f};
        normal = multiplier * computeFaceNormal(u, v, w);
        // dupliacte vertices

        // extrude vertex
        vertexPos[V] = vertexPos[V] + normal;

        // update variables
        u = v;
        v = w;

        x = halfEdge[h][NEXT];
        halfEdge[h][NEXT] = H;
        h = x;

    } while (h != start);
}

void MeshHalfEdge::addFace(const Face &face, IntList &prev)
{
    int n = face.size();
    assert(n > 0);
    Key key;
    // HalfEdge *h;

    int faceIdx = FaceHalfEdge.size();

    int m = halfEdge.size(); // Get starting index for new half-edges

    halfEdge.resize(m + n, {-1, -1, -1, -1}); // Resize for new half-edges
    prev.resize(m + n, -1);                   // Resize for new half-edges
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
        // std::cerr << "setting half edge for vertex" << vertex << "to " << h << std::endl;
        // vertexHalfEdge[vertex] = h;
        if (vertexHalfEdge[vertex] == -1)
        {
            vertexHalfEdge[vertex] = h;
        }
        prev[h] = vertex;

        key = enkey(nextVertex, vertex);

        if (halfedgeMap.find(key) != halfedgeMap.end())
        {
            print("Found key ", key, "For vertices", nextVertex, vertex);
            // print("")
            pair = halfedgeMap[key];
            halfEdge[h][PAIR] = pair;
            halfEdge[pair][PAIR] = h;
        }
        else
        {
            key = enkey(vertex, nextVertex);

            print("Added key ", key, "For vertices", vertex, nextVertex);
            halfedgeMap[key] = h;
        }

        // halfEdge[h][PREV] = (i - 1 + n) % n + m;

        halfEdge[h][NEXT] = (i + 1) % n + m;

        halfEdge[h][LEFT] = faceIdx;
        halfEdge[h][HEAD] = nextVertex;
    }
    int flag = 1;
    for (int i = 0; i < n; i++)
    {
        int v = face[i];
        // if (v == 13)
        // {
        //     std::cerr << "Assigned half edges to " << v << "as" << vertexHalfEdge[v] << std::endl;
        // }
        if (vertexHalfEdge[v] == -1)
        {
            flag = 0;
            std::cerr << "ERROR HALF EDGE NOT SET FOR " << v << std::endl;
        }
        // if (vertexHalfEdge[v] == -1) {
        //     std::cerr << "ERROR: Vertex " << v << " is not connected to any half-edge!" << std::endl;
        // }
    }
    if (flag)
    {
        // std::cerr << "All these vertices have half edges :";
        // for (auto i : face)
        //     std::cerr << i << " , ";
        // std::cerr << std::endl;
    }
};

void MeshHalfEdge::handleBoundaryVertices(IntList &prev)
{
    // Define : boudary half edge : pair == -1
    // Assumption : a vertex has either no boundary half edges, or exactly one outgoign and one incoming boundary half edge .

    // using HalfEdge = std::array<int, 4>;
    // vertexhalfedge -> outgoing half edges
    int u, v;
    int n = halfEdge.size();
    for (int h = 0; h < n; h++)
    {
        if (halfEdge[h][PAIR] != -1)
            continue;
        u = prev[h];
        v = halfEdge[h][HEAD];
        halfEdge[h][PAIR] = halfEdge.size();
        vertexHalfEdge[v] = halfEdge.size();
        HalfEdge dummy;
        dummy[HEAD] = u;
        dummy[PAIR] = h;
        dummy[NEXT] = -1;
        dummy[LEFT] = -1;
        halfEdge.push_back(dummy);
    }
}

void MeshHalfEdge::buildHalfEdgeStructure(VecList &faces)
{
    // Assumption : either manifold or manifold with boundary
    // Hence : No half edges has NEXT == -1
    initArray();
    // IntMap boundaryEdgeTailMap; // Maps boundary edge h -> tail vertex u
    IntList prev;

    for (Face &face : faces)
    {
        addFace(face, prev);
    }
    handleBoundaryVertices(prev);
    if (vertexNormal.size() < vertexPos.size())
    {
        vertexNormal.resize(vertexPos.size(), Vec3(0.0f, 0.0f, 0.0f)); // Extend with zero normals
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
    IntList prev;
    for (auto &tri : triangles)
    {
        addFace(triangle_to_face(tri), prev);
    }
    handleBoundaryVertices(prev);
    if (vertexNormal.size() < vertexPos.size())
    {
        vertexNormal.resize(vertexPos.size(), Vec3(0.0f, 0.0f, 0.0f)); // Extend with zero normals
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
    print("Triangulate Face called");
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
            print("Added edge 1 : ", x, "->", z);
        }
        else
        {
            // edges.push_back({y, z});
            edges.push_back(Edge(y, z));
            print("Added edge 2 : ", y, "->", z);

            IVec3 tri(x, y, z);
            // IVec3(x, y, z)
            triangleVertices.push_back(tri);
            // computeTriangleNormal(tri);

            // faceNormals
        }

        y = z;
        h = halfEdge[h[NEXT]];
        print("x,y,z", x, y, z);
    }

    edges.push_back(Edge(z, x));
    print("Added edge 2 : ", z, "->", x);
}

// void MeshHalfEdge::computeTriangleNormal(IVec3 &tri)
// {

//     glm::vec3 edge1 = vertexPos[tri.z] - vertexPos[tri.y];
//     glm::vec3 edge2 = vertexPos[tri.y] - vertexPos[tri.x];

//     glm::vec3 faceNormal = glm::normalize(glm::cross(edge1, edge2));
//     faceNormals.push_back(faceNormal);
// }

Vec3 MeshHalfEdge::computeFaceNormal(int x, int y, int z)
{

    glm::vec3 edge1 = vertexPos[z] - vertexPos[y];
    glm::vec3 edge2 = vertexPos[y] - vertexPos[x];

    glm::vec3 faceNormal = -glm::normalize(glm::cross(edge1, edge2));
    // faceNormals.push_back(faceNormal);
    return faceNormal;
}
// void getAdjacentFaces()
// {
// }

// void computeVertexNormals_helper()
// {
//     // FaceList
//     IntList faces = getAdjacentFaces();
//     Vec3 normal;
//     for (idx in faces)
//     {
//         normal = normal + faceNormals[idx];
//     }
//     normal = normal / faces.size();
// }

void MeshHalfEdge::sanity_check()
{
    std::cout << "Running Half-Edge Sanity Check..." << std::endl;

    for (int h = 0; h < halfEdge.size(); h++)
    {
        int next = halfEdge[h][NEXT];
        int pair = halfEdge[h][PAIR];
        int head = halfEdge[h][HEAD];

        // Check valid NEXT pointer
        if (next == -1)
        {
            std::cerr << "ERROR: Half-edge " << h << " has an invalid NEXT pointer!" << std::endl;
        }

        if (pair == -1)
        {
            std::cerr << "ERROR: Half-edge " << h << " has an invalid PAIR  pointer!" << std::endl;
        }
        else if (halfEdge[pair][PAIR] != h)
        // if (pair != -1 && halfEdge[pair][PAIR] != h)
        {
            std::cerr << "ERROR: Half-edge " << h << " and its pair " << pair << " (points to " << halfEdge[pair][PAIR] << ") are not mutual pairs!" << std::endl;
        }

        // Check valid HEAD
        if (head < 0 || head >= vertexPos.size())
        {
            std::cerr << "ERROR: Half-edge " << h << " points to an invalid vertex : " << head << std::endl;
        }
    }

    for (int v = 0; v < vertexPos.size(); v++)
    {
        if (vertexHalfEdge[v] == -1)
        {
            std::cerr << "WARNING: Vertex " << v << " is isolated (not connected to any half-edge)!" << std::endl;
        }
    }

    std::cout << "Checking vertex connectivity..." << std::endl;

    for (int v = 0; v < vertexPos.size(); v++)
    {
        int h = vertexHalfEdge[v]; // Get outgoing half-edge from vertex
        if (h == -1)
        {
            std::cerr << "WARNING: Vertex " << v << " is isolated!" << std::endl;
            continue;
        }

        int start = h;
        int count = 0;

        // std::cout << "Neighbors of vertex " << v << ": ";

        // do
        // {
        //     if (h == -1)
        //     {
        //         std::cerr << "ERROR: Null half-edge encountered at vertex " << v << "!" << std::endl;
        //         break;
        //     }

        //     int neighbor = halfEdge[h][HEAD];
        //     std::cout << neighbor << " ";

        //     h = halfEdge[h][PAIR]; // Move to opposite half-edge

        //     if (h == -1)
        //     {
        //         std::cerr << "ERROR: Boundary edge encountered at vertex " << v << "!" << std::endl;
        //         break;
        //     }

        //     h = halfEdge[h][NEXT]; // Move to next outgoing edge

        //     count++;

        //     if (count > 20) // Detect infinite loop
        //     {
        //         std::cerr << "ERROR: Possible infinite loop at vertex " << v << "!" << std::endl;
        //         break;
        //     }

        // } while (h != start); // Should complete the loop

        std::cout << std::endl;
    }

    std::cout << "Vertex neighbor check complete!" << std::endl;

    std::cerr << "Sanity Check Complete!" << std::endl;
}

void MeshHalfEdge::computeVertexNormals()
{
    int h, f, count;
    // Vec3 normal;

    int start;
    int u, w;
    for (int v = 0; v < vertexPos.size(); v++)
    {
        // std::cerr << "Recomputing normals for vertex" << v << std::endl;

        Vec3 normal = {0.0f, 0.0f, 0.0f};
        h = vertexHalfEdge[v];

        // h is outgoing half edge from vertex v

        start = h;
        count = 0;
        do
        {
            // assert(halfEdge[h]!=)
            // printerr("-------------------");
            u = halfEdge[h][HEAD];
            // printerr("CURRENT HALF EDGE : ", h, "u = ", u);
            h = halfEdge[h][PAIR];
            // printerr("CURRENT HALF EDGE after pair : ", h);
            // printerr("CURRENT HALF EDGE : ", h);

            // assert(h != -1);

            assert(halfEdge[h][HEAD] == v);

            // if (count > 8)
            // {
            //     print("INFINITE LOOP ERROR");
            //     printerr("INFINITE LOOP ERROR");
            //     break;
            // }
            // count++;
            if (halfEdge[h][NEXT] == -1)
            {
                // Must be a boundary edge, so LEFT and NEXT both must be 0
                assert(halfEdge[h][LEFT] == -1);
                // printerr("This is a boundary edge:", h);
                break;
            }

            h = halfEdge[h][NEXT];

            w = halfEdge[h][HEAD];
            // printerr("CURRENT HALF EDGE after next : ", h, "W = ", w);

            normal = normal + computeFaceNormal(u, v, w);
            // print(u, "->", v, "->", w);
        } while (h != start);
        vertexNormal[v] = glm::normalize(normal);
        // vertexNormal[v] = normal / (float)count;
    }
}

void zeroVec(Vec3 &vec)
{
    vec.x = 0.0f;
    vec.y = 0.0f;
    vec.z = 0.0f;
}

void MeshHalfEdge::smoothen_step(float λ)
{
    int start;
    int h;

    int v1;
    float count;
    for (int v = 0; v < vertexPos.size(); v++)
    {
        // zeroVec(sum);
        if (false)
            print("Smoothen vertex", v);
        Vec3 sum = Vec3(0.0f);

        count = 0;
        start = vertexHalfEdge[v];
        h = start;
        do
        {
            count++;
            if (false)
                print(count, "number of neighbours");
            v1 = halfEdge[h][HEAD];
            sum = sum + vertexPos[v1];
            h = halfEdge[h][PAIR];
            if (halfEdge[h][NEXT] == -1)
            {
                // Must be a boundary edge, so LEFT and NEXT both must be 0
                assert(halfEdge[h][LEFT] == -1);
                // printerr("This is a boundary edge:", h);
                break;
            }
            h = halfEdge[h][NEXT];

        } while (h != start);

        if (count == 0.0f)
            continue;
        if (false)
        {
            print("count is", count);
            print("sum is :", sum.x, " ", sum.y, " ", sum.z);
        }
        sum = sum / count;
        sum = sum - vertexPos[v];
        // print("updated by )
        vertexPos[v] = vertexPos[v] + sum * λ;
    }
}
void MeshHalfEdge::smoothen_taubin(float λ, float u, int iterations)
{
    // taubin
    // traverse all vertices adjacent to a given vertex
    // Vec3 sum(0.0f, 0.0f, 0.0f);
    assert(λ >= 0);
    assert(λ <= 1);
    assert(u >= -11);
    assert(u <= 0);

    for (int iter = 0; iter < iterations; iter++)
    {
        smoothen_step(λ);
        smoothen_step(u);
    }
}

void MeshHalfEdge::smoothen_laplacian(float λ, int iterations)
{
    // taubin
    // traverse all vertices adjacent to a given vertex
    // Vec3 sum(0.0f, 0.0f, 0.0f);
    smoothen_taubin(λ, 0.0f, iterations);
}

// void MeshHalfEdge::extractEdgesFromFaces();

// // part2

void saveOBJ(const std::string &filename, const std::vector<Vec3> &vertices, const std::vector<Vec3> &normals, const std::vector<std::vector<int>> &faces)
{
    std::ofstream file(filename);
    if (!file)
    {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    for (const auto &v : vertices)
        file << "v " << v.x << " " << v.y << " " << v.z << "\n";

    for (const auto &n : normals)
        file << "vn " << n.x << " " << n.y << " " << n.z << "\n";

    for (const auto &f : faces)
    {
        file << "f";
        for (int idx : f)
            file << " " << (idx + 1) << "//" << (idx + 1); // OBJ indices are 1-based
        file << "\n";
    }
}

int sign(float x)
{
    if (x >= 0)
        return 1;
    else
        return -1;
}

void generateGrid_helper(int m, int n, float a, int axis,
                         std::vector<Vec3> &vertices,
                         std::vector<Vec3> &normals,
                         std::vector<std::vector<int>> &faces)
{
    int startIdx = vertices.size(); // Track starting index for face indexing
    float k = sign(a);

    Vec3 normal;
    if (axis == 0)
        normal = {k, 0.0f, 0.0f}; // Plane at x = a
    else if (axis == 1)
        normal = {0.0f, k, 0.0f}; // Plane at y = a
    else if (axis == 2)
        normal = {0.0f, 0.0f, k}; // Plane at z = a
    else
        return; // Invalid axis, do nothing

    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            Vec3 vertex;
            if (axis == 0)
            { // x = a
                vertex = {a, (float)i / m - 0.5f, (float)j / n - 0.5f};
            }
            else if (axis == 1)
            { // y = a
                vertex = {(float)j / n - 0.5f, a, (float)i / m - 0.5f};
            }
            else if (axis == 2)
            { // z = a
                vertex = {(float)i / m - 0.5f, (float)j / n - 0.5f, a};
                // vertex = {(float)j / n - 0.5f, (float)i / m - 0.5f, a};
            }
            vertices.push_back(vertex);
            normals.push_back(normal);
        }
    }

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int idx1 = startIdx + i * (n + 1) + j;
            int idx2 = idx1 + 1;
            int idx3 = idx1 + (n + 1);
            int idx4 = idx3 + 1;
            if (((k < 0) && (axis == 2)) || ((k > 0) && (axis != 2)))
                faces.push_back({idx1, idx2, idx4, idx3});
            else
                faces.push_back({idx3, idx4, idx2, idx1});
        }
    }
}

void generateCustomGrid(int m, int n, float a, int axis, const std::string &filename)
{
    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    std::vector<std::vector<int>> faces;

    generateGrid_helper(m, n, a, axis,
                        vertices,
                        normals,
                        faces);

    saveOBJ(filename, vertices, normals, faces);
}

void oldmakecube1(int m, int n, int o, const std::string &filename)
{
    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    std::vector<std::vector<int>> faces;

    int axis = 3;

    float a = 1.5;
    do
    {
        a--;

        generateGrid_helper(n, o, a, 0,
                            vertices,
                            normals,
                            faces);

        generateGrid_helper(o, m, a, 1,
                            vertices,
                            normals,
                            faces);

        generateGrid_helper(m, n, a, 2,
                            vertices,
                            normals,
                            faces);
    } while (a >= 0);

    saveOBJ(filename, vertices, normals, faces);
}

void oldmakecube2(int m, int n, int o, const std::string &filename)
{
    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    std::vector<std::vector<int>> faces;

    // Precompute unique vertices
    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            for (int k = 0; k <= o; k++)
            {
                Vec3 vertex = {(float)i / m - 0.5f, (float)j / n - 0.5f, (float)k / o - 0.5f};
                vertices.push_back(vertex);
            }
        }
    }

    auto getIndex = [&](int i, int j, int k) -> int
    {
        return i * (n + 1) * (o + 1) + j * (o + 1) + k;
    };

    auto addFace = [&](int i1, int j1, int k1,
                       int i2, int j2, int k2,
                       int i3, int j3, int k3,
                       int i4, int j4, int k4)
    {
        faces.push_back({getIndex(i1, j1, k1), getIndex(i2, j2, k2),
                         getIndex(i3, j3, k3), getIndex(i4, j4, k4)});
    };

    // ** Z = 0 (Back Face, Normal (0, 0, -1)) **
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            addFace(i, j, 0, i, j + 1, 0, i + 1, j + 1, 0, i + 1, j, 0); // CCW

    // ** Z = o (Front Face, Normal (0, 0, +1)) **
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            addFace(i, j, o, i + 1, j, o, i + 1, j + 1, o, i, j + 1, o); // CCW

    // ** X = 0 (Left Face, Normal (-1, 0, 0)) **
    for (int j = 0; j < n; j++)
        for (int k = 0; k < o; k++)
            addFace(0, j, k, 0, j, k + 1, 0, j + 1, k + 1, 0, j + 1, k); // CCW

    // ** X = m (Right Face, Normal (+1, 0, 0)) **
    for (int j = 0; j < n; j++)
        for (int k = 0; k < o; k++)
            addFace(m, j, k, m, j + 1, k, m, j + 1, k + 1, m, j, k + 1); // CCW

    // ** Y = 0 (Bottom Face, Normal (0, -1, 0)) **
    for (int i = 0; i < m; i++)
        for (int k = 0; k < o; k++)
            addFace(i, 0, k, i + 1, 0, k, i + 1, 0, k + 1, i, 0, k + 1); // CCW

    // ** Y = n (Top Face, Normal (0, +1, 0)) **
    for (int i = 0; i < m; i++)
        for (int k = 0; k < o; k++)
            addFace(i, n, k, i, n, k + 1, i + 1, n, k + 1, i + 1, n, k); // CCW

    // Precompute unique normals for each vertex
    normals.resize(vertices.size(), {0, 0, 0});

    auto setNormal = [&](int i, int j, int k, Vec3 normal)
    {
        int idx = getIndex(i, j, k);
        if (normals[idx] == Vec3{0, 0, 0}) // Assign only if not set
            normals[idx] = normal;
    };

    // Set normals for each face's vertices
    for (int i = 0; i <= m; i++)
        for (int j = 0; j <= n; j++)
            for (int k = 0; k <= o; k++)
            {
                if (i == 0)
                    setNormal(i, j, k, {-1, 0, 0});
                if (i == m)
                    setNormal(i, j, k, {1, 0, 0});
                if (j == 0)
                    setNormal(i, j, k, {0, -1, 0});
                if (j == n)
                    setNormal(i, j, k, {0, 1, 0});
                if (k == 0)
                    setNormal(i, j, k, {0, 0, -1});
                if (k == o)
                    setNormal(i, j, k, {0, 0, 1});
            }
    saveOBJ(filename, vertices, normals, faces);
}
void generateCube_3(int m, int n, int o, const std::string &filename)
{
    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    std::vector<std::vector<int>> faces;

    // Precompute unique vertices
    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            for (int k = 0; k <= o; k++)
            {
                Vec3 vertex = {(float)i / m - 0.5f, (float)j / n - 0.5f, (float)k / o - 0.5f};
                vertices.push_back(vertex);
            }
        }
    }

    auto getIndex = [&](int i, int j, int k) -> int
    {
        assert(i <= n);
        auto x = i * (n + 1) * (o + 1) + j * (o + 1) + k;
        if (x == 13 | x == 12)
        {
            std::cerr << "Accessed x : " << x << "i,j,k" << i << " " << j << " " << k << std::endl;
        }
        return x;
        // return i * (n + 1) * (o + 1) + j * (o + 1) + k;
    };

    auto addFace = [&](int i1, int j1, int k1,
                       int i2, int j2, int k2,
                       int i3, int j3, int k3,
                       int i4, int j4, int k4)
    {
        faces.push_back({getIndex(i1, j1, k1), getIndex(i2, j2, k2), getIndex(i3, j3, k3), getIndex(i4, j4, k4)});
    };

    // ** Z = 0 (Back Face, Normal (0, 0, -1)) **
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            addFace(i, j, 0, i, j + 1, 0, i + 1, j + 1, 0, i + 1, j, 0); // CCW

    // ** Z = o (Front Face, Normal (0, 0, +1)) **
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            addFace(i, j, o, i + 1, j, o, i + 1, j + 1, o, i, j + 1, o); // CCW

    // ** X = 0 (Left Face, Normal (-1, 0, 0)) **
    for (int j = 0; j < n; j++)
        for (int k = 0; k < o; k++)
            addFace(0, j, k, 0, j, k + 1, 0, j + 1, k + 1, 0, j + 1, k); // CCW

    // ** X = m (Right Face, Normal (+1, 0, 0)) **
    for (int j = 0; j < n; j++)
        for (int k = 0; k < o; k++)
            addFace(m, j, k, m, j + 1, k, m, j + 1, k + 1, m, j, k + 1); // CCW

    // ** Y = 0 (Bottom Face, Normal (0, -1, 0)) **
    for (int i = 0; i < m; i++)
        for (int k = 0; k < o; k++)
            addFace(i, 0, k, i + 1, 0, k, i + 1, 0, k + 1, i, 0, k + 1); // CCW

    // ** Y = n (Top Face, Normal (0, +1, 0)) **
    for (int i = 0; i < m; i++)
        for (int k = 0; k < o; k++)
            addFace(i, n, k, i, n, k + 1, i + 1, n, k + 1, i + 1, n, k); // CCW

    // Precompute unique normals for each vertex
    normals.resize(vertices.size(), {0, 0, 0});

    auto setNormal = [&](int i, int j, int k, Vec3 normal)
    {
        int idx = getIndex(i, j, k);
        if (normals[idx] == Vec3{0, 0, 0}) // Assign only if not set
            normals[idx] = normal;
    };

    // Set normals for each face's vertices
    for (int i = 0; i <= m; i++)
        for (int j = 0; j <= n; j++)
            for (int k = 0; k <= o; k++)
            {
                if (i == 0)
                    setNormal(i, j, k, {-1, 0, 0});
                if (i == m)
                    setNormal(i, j, k, {1, 0, 0});
                if (j == 0)
                    setNormal(i, j, k, {0, -1, 0});
                if (j == n)
                    setNormal(i, j, k, {0, 1, 0});
                if (k == 0)
                    setNormal(i, j, k, {0, 0, -1});
                if (k == o)
                    setNormal(i, j, k, {0, 0, 1});
            }
    saveOBJ(filename, vertices, normals, faces);
}
void generateCube(int m, int n, int o, const std::string &filename)
{
    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    std::vector<std::vector<int>> faces;

    int total = (m + 1) * (n + 1) * (o + 1);
    // std::map<std::tuple<int, int, int>, int> vertexMap; // Map to store unique surface vertices
    IntList verticeMap(total, -1);
    // ** Insert only SURFACE vertices **

    // auto insertVertex = [&](int i, int j, int k) -> int
    // {
    //     std::tuple<int, int, int> key = {i, j, k};
    //     if (vertexMap.count(key))
    //         return vertexMap[key]; // Already exists

    //     Vec3 vertex = {(float)i / m - 0.5f, (float)j / n - 0.5f, (float)k / o - 0.5f};
    //     int index = vertices.size();
    //     vertices.push_back(vertex);
    //     vertexMap[key] = index;
    //     return index;
    // };

    auto setNormalForFace = [&](int i, int j, int k)
    {
        if (i == 0)
            normals.push_back({-1, 0, 0});
        else if (i == m)
            normals.push_back({1, 0, 0});
        else if (j == 0)
            normals.push_back({0, -1, 0});
        else if (j == n)
            normals.push_back({0, 1, 0});
        else if (k == 0)
            normals.push_back({0, 0, -1});
        else if (k == o)
            normals.push_back({0, 0, 1});
    };

    auto getIndex = [&](int i, int j, int k) -> int
    {
        // assert(i <= n);
        auto x = i * (n + 1) * (o + 1) + j * (o + 1) + k;
        if (verticeMap[x] == -1)
        {
            verticeMap[x] = vertices.size();
            Vec3 vertex = {(float)i / m - 0.5f, (float)j / n - 0.5f, (float)k / o - 0.5f};
            vertices.push_back(vertex);
            setNormalForFace(i, j, k);
        }

        return verticeMap[x];
    };

    // auto setNormal = [&](int idx, Vec3 normal)
    // {

    //         // normals[idx] = normal;
    // };

    auto addFace = [&](int i1, int j1, int k1, int i2, int j2, int k2, int i3, int j3, int k3, int i4, int j4, int k4)
    {
        faces.push_back({getIndex(i1, j1, k1), getIndex(i2, j2, k2), getIndex(i3, j3, k3), getIndex(i4, j4, k4)});
        // setNormalForFace(i1, j1, k1);
        // setNormalForFace(i2, j2, k2);
        // setNormalForFace(i3, j3, k3);
        // setNormalForFace(i4, j4, k4);
    };

    // ** Z = 0 (Back Face, Normal (0, 0, -1)) **
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            addFace(i, j, 0, i, j + 1, 0, i + 1, j + 1, 0, i + 1, j, 0); // CCW

    // ** Z = o (Front Face, Normal (0, 0, +1)) **
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            addFace(i, j, o, i + 1, j, o, i + 1, j + 1, o, i, j + 1, o); // CCW

    // ** X = 0 (Left Face, Normal (-1, 0, 0)) **
    for (int j = 0; j < n; j++)
        for (int k = 0; k < o; k++)
            addFace(0, j, k, 0, j, k + 1, 0, j + 1, k + 1, 0, j + 1, k); // CCW

    // ** X = m (Right Face, Normal (+1, 0, 0)) **
    for (int j = 0; j < n; j++)
        for (int k = 0; k < o; k++)
            addFace(m, j, k, m, j + 1, k, m, j + 1, k + 1, m, j, k + 1); // CCW

    // ** Y = 0 (Bottom Face, Normal (0, -1, 0)) **
    for (int i = 0; i < m; i++)
        for (int k = 0; k < o; k++)
            addFace(i, 0, k, i + 1, 0, k, i + 1, 0, k + 1, i, 0, k + 1); // CCW

    // ** Y = n (Top Face, Normal (0, +1, 0)) **
    for (int i = 0; i < m; i++)
        for (int k = 0; k < o; k++)
            addFace(i, n, k, i, n, k + 1, i + 1, n, k + 1, i + 1, n, k); // CCW

    // Precompute unique normals for each vertex

    // normals.resize(()).size(), {0, 0, 0});

    // // Set normals for each face's vertices
    // for (int i = 0; i <= m; i++)
    //     for (int j = 0; j <= n; j++)
    //         for (int k = 0; k <= o; k++)
    //         {
    //             if (i == 0)
    //                 setNormal(i, j, k, {-1, 0, 0});
    //             if (i == m)
    //                 setNormal(i, j, k, {1, 0, 0});
    //             if (j == 0)
    //                 setNormal(i, j, k, {0, -1, 0});
    //             if (j == n)
    //                 setNormal(i, j, k, {0, 1, 0});
    //             if (k == 0)
    //                 setNormal(i, j, k, {0, 0, -1});
    //             if (k == o)
    //                 setNormal(i, j, k, {0, 0, 1});
    //         }
    // std::cerr << "Vertices"<<filename
    printerr("Actual Vertices : ", vertices.size(), " | Normals : ", normals.size(), " | Faces : ", faces.size());
    int exp = 8 + 4 * (m - 1 + n - 1 + o - 1) + 2 * ((m - 1) * (n - 1) + (n - 1) * (o - 1) + (o - 1) * (m - 1));
    printerr("Expect Vertices : ", exp, " | Normals : ", exp, " | Faces : ", 2 * (m * n + n * o + o * m));
    saveOBJ(filename, vertices, normals, faces);
}

void generateGrid(int m, int n, const std::string &filename)
{
    generateCustomGrid(m, n, 0.0f, 2,
                       filename);
}

// // void generateGrid(int m, int n, const std::string &filename)
// {
//     std::vector<Vec3> vertices;
//     std::vector<Vec3> normals;
//     std::vector<std::vector<int>> faces;

//     for (int i = 0; i <= m; i++)
//     {
//         for (int j = 0; j <= n; j++)
//         {
//             vertices.push_back({(float)j / n, (float)i / m, 0.0f});
//             normals.push_back({0.0f, 0.0f, -1.0f});
//         }
//     }

//     for (int i = 0; i < m; i++)
//     {
//         for (int j = 0; j < n; j++)
//         {
//             int idx1 = i * (n + 1) + j;
//             int idx2 = idx1 + 1;
//             int idx3 = idx1 + (n + 1);
//             int idx4 = idx3 + 1;
//             faces.push_back({idx1, idx2, idx4, idx3});
//         }
//     }

//     saveOBJ(filename, vertices, normals, faces);
// }

// void pushSphereFace(FaceList &faces, )

void generateSphere(int m, int n, const std::string &filename)

{
    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    std::vector<std::vector<int>> faces;

    // auto addFace = [&](int a, int b, int c)
    // {
    //     faces.push_back({a, b, c});
    // };

    auto addFace = [&](std::initializer_list<int> indices)
    {
        std::vector<int> face(indices);
        std::reverse(face.begin(), face.end()); // Reverse before inserting
        faces.push_back(face);
    };

    vertices.push_back({0, 0, 1}); // Top pole
    normals.push_back({0, 0, 1});

    for (int i = 1; i < n; i++)
    {
        float phi = M_PI * i / n;
        for (int j = 0; j < m; j++)
        {
            float theta = 2 * M_PI * j / m;
            float x = sin(phi) * cos(theta);
            float y = sin(phi) * sin(theta);
            float z = cos(phi);
            vertices.push_back({x, y, z});
            normals.push_back({x, y, z});
        }
    }

    vertices.push_back({0, 0, -1}); // Bottom pole
    normals.push_back({0, 0, -1});

    for (int j = 0; j < m; j++)
    {
        // faces.push_back({(j + 1) % m + 1, j + 1, 0});
        // faces.push_back({(j + 1) % m + 1, j + 1, 0});
        addFace({(j + 1) % m + 1, j + 1, 0});
        // addFace((j + 1) % m + 1, j + 1, 0);
        // faces.push_back({0, j + 1, (j + 1) % m + 1});
    }

    for (int i = 0; i < n - 2; i++)
    {
        for (int j = 0; j < m; j++)
        {
            int idx1 = i * m + j + 1;
            int idx2 = idx1 + 1;
            int idx3 = idx1 + m;
            int idx4 = idx3 + 1;
            if ((j + 1) % m == 0)
            {
                idx2 -= m;
                idx4 -= m;
            }
            // faces.push_back({idx1, idx2, idx4, idx3});
            // addFace(idx1, idx2, idx4, idx3);
            addFace({idx1, idx2, idx4, idx3});
            // faces.push_back({idx4, idx3, idx2, idx1});
        }
    }

    int last = vertices.size() - 1;
    for (int j = 0; j < m; j++)
    {
        // addFace(last, last - m + j, last - m + (j + 1) % m);
        addFace({last, last - m + j, last - m + (j + 1) % m});
        // faces.push_back({last, last - m + j, last - m + (j + 1) % m});
        // faces.push_back({last - m + (j + 1) % m, last - m + j, last});
    }

    saveOBJ(filename, vertices, normals, faces);
}

// void generateCube(int m, int n, int o, const std::string &filename)
// {
//     std::vector<Vec3> vertices;
//     std::vector<Vec3> normals;
//     std::vector<std::vector<int>> faces;

//     for (int i = 0; i <= m; i++)
//     {
//         float x = -0.5f + (float)i / m;
//         for (int j = 0; j <= n; j++)
//         {
//             float y = -0.5f + (float)j / n;
//             for (int k = 0; k <= o; k++)
//             {
//                 float z = -0.5f + (float)k / o;
//                 vertices.push_back({x, y, z});
//             }
//         }
//     }

//     normals = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};

//     // Generate faces for the six sides of the cube (both front and back)
//     for (int i = 0; i < m; i++)
//     {
//         for (int j = 0; j < n; j++)
//         {
//             int idx1 = i * (n + 1) * (o + 1) + j * (o + 1);
//             int idx2 = idx1 + (o + 1);
//             int idx3 = idx1 + 1;
//             int idx4 = idx2 + 1;
//             faces.push_back({idx1, idx3, idx4, idx2});
//             faces.push_back({idx2, idx4, idx3, idx1});
//         }
//     }
//     for (int j = 0; j < n; j++)
//     {
//         for (int k = 0; k < o; k++)
//         {
//             int idx1 = j * (o + 1) + k;
//             int idx2 = idx1 + 1;
//             int idx3 = idx1 + (n + 1) * (o + 1);
//             int idx4 = idx3 + 1;
//             faces.push_back({idx1, idx2, idx4, idx3});
//             faces.push_back({idx3, idx4, idx2, idx1});
//         }
//     }
//     for (int i = 0; i < m; i++)
//     {
//         for (int k = 0; k < o; k++)
//         {
//             int idx1 = i * (n + 1) * (o + 1) + k;
//             int idx2 = idx1 + 1;
//             int idx3 = idx1 + (n + 1) * (o + 1);
//             int idx4 = idx3 + 1;
//             faces.push_back({idx1, idx2, idx4, idx3});
//             faces.push_back({idx3, idx4, idx2, idx1});
//         }
//     }

//     saveOBJ(filename, vertices, normals, faces);
//     // std::ofstream file(filename);
//     // if (!file) {
//     //     std::cerr << "Failed to open file: " << filename << "\n";
//     //     return;
//     // }

//     // for (const auto& v : vertices)
//     //     file << "v " << v.x << " " << v.y << " " << v.z << "\n";

//     // for (const auto& n : normals)
//     //     file << "vn " << n.x << " " << n.y << " " << n.z << "\n";

//     // for (const auto& f : faces) {
//     //     file << "f";
//     //     for (int idx : f) file << " " << (idx + 1) << "//" << ((idx % 6) + 1);
//     //     file << "\n";
//     // }
// }

// void MeshHalfEdge::generateGridMesh(int m, int n);
// void MeshHalfEdge::generateSphereMesh(int m, int n);
// void MeshHalfEdge::enerateCubeMesh(int m, int n, int o);

// // part3
// void MeshHalfEdge::loadOBJ(std::string filename);
// void MeshHalfEdge::saveOBJ(std::string filename);

// // part4
// void MeshHalfEdge::computeFaceNormals();
// void MeshHalfEdge::computeVertexNormals();