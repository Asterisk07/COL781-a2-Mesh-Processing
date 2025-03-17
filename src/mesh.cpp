#include "mesh.hpp" // ✅ Always include the corresponding header
#include <glm/glm.hpp>
#include <random>
#include <string>

// Helper functions

int sign(float x)
{
    if (x >= 0)
        return 1;
    else
        return -1;
}

Key enkey(int x, int y)
{
    return ((Key(x) << 32) | y);
}

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

// part1

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
        // std::cout << "vertex is " << vertex << std::endl;
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
            ////print("Found key ", key, "For vertices", nextVertex, vertex);
            pair = halfedgeMap[key];
            halfEdge[h][PAIR] = pair;
            halfEdge[pair][PAIR] = h;
        }
        else
        {
            key = enkey(vertex, nextVertex);

            ////print("Added key ", key, "For vertices", vertex, nextVertex);
            halfedgeMap[key] = h;
        }

        // halfEdge[h][PREV] = (i - 1 + n) % n + m;

        halfEdge[h][NEXT] = (i + 1) % n + m;

        halfEdge[h][LEFT] = faceIdx;
        halfEdge[h][HEAD] = nextVertex;
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
    vertexHalfEdge.resize(vertexPos.size(), -1);
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

void MeshHalfEdge::triangulateFace(int faceIdx, IVec3List &triangleVertices, EdgeList &edges)
{
    // triangulates face and adds triangle to triangleVertices
    // int h = FaceHalfEdge[faceIdx];
    // HalfEdge h =
    // print("Triangulate Face called");
    HalfEdge h;

    h = halfEdge[FaceHalfEdge[faceIdx]];

    int x = -1, y = -1, z = -1;
    x = h[HEAD];
    // print("currently h[HEAD] = ", h[HEAD]);

    h = halfEdge[h[NEXT]];
    // print("x,y,z", x, y, z);

    while (h[HEAD] != x)
    {
        z = h[HEAD];
        if (y == -1)
        {
            if (x < z)
                edges.push_back(Edge(x, z));
            // printerr("Added edge 1 : ", x, "->", z);
        }
        else
        {
            // edges.push_back({y, z});
            if (y < z)
                edges.push_back(Edge(y, z));
            // printerr("Added edge 2 : ", y, "->", z);

            IVec3 tri(x, y, z);
            // IVec3(x, y, z)
            triangleVertices.push_back(tri);

            // faceNormals
        }

        y = z;
        h = halfEdge[h[NEXT]];
        // print("x,y,z", x, y, z);
    }

    if (z < x)
        edges.push_back(Edge(z, x));
    // printerr("Added edge 3 : ", z, "->", x);
}

void MeshHalfEdge::triangulateMesh(IVec3List &triangleVertices, EdgeList &edges)
{
    int numFaces = FaceHalfEdge.size();
    for (int faceidx = 0; faceidx < numFaces; faceidx++)
    {
        triangulateFace(faceidx, triangleVertices, edges);
    }
}

void MeshHalfEdge::sanity_check()
{
    // std::cout << "Running Half-Edge Sanity Check..." << std::endl;

    for (int h = 0; h < halfEdge.size(); h++)
    {
        int next = halfEdge[h][NEXT];
        int pair = halfEdge[h][PAIR];
        int head = halfEdge[h][HEAD];
        int face = halfEdge[h][LEFT];

        // Check valid NEXT pointer
        if (next == -1)
        {
            assert(face == -1);
            // std::cerr << "BOUNDARY HALF : Half-edge " << h << " has an invalid NEXT pointer!" << std::endl;
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

    // std::cout << "Checking vertex connectivity..." << std::endl;

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

        // std::cout << std::endl;
    }

    // std::cout << "Vertex neighbor check complete!" << std::endl;

    std::cerr << "Sanity Check Complete!" << std::endl;
}

void MeshHalfEdge::debugInfo(IVec3List &triangleVertices, EdgeList &edges)
{

    std::cerr << " Mesh info :  Vertices " << vertexPos.size() << " | Half Edges " << halfEdge.size() << " | Faces " << FaceHalfEdge.size() << "\n Rendering info : triangles " << triangleVertices.size() << "  Edges " << edges.size() << std::endl;
}

// part 2
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

    for (int j = 0; j <= n; j++)
    {
        for (int i = 0; i <= m; i++)
        {
            Vec3 vertex;
            if (axis == 0)
            { // x = a
                vertex = {a, (float)j / n - 0.5f, (float)i / m - 0.5f};
            }
            else if (axis == 1)
            { // y = a
                vertex = {(float)i / m - 0.5f, a, (float)j / n - 0.5f};
            }
            else if (axis == 2)
            { // z = a
                vertex = {(float)j / n - 0.5f, (float)i / m - 0.5f, a};
                // vertex = {(float)i / m - 0.5f, (float)j / n - 0.5f, a};
            }
            vertices.push_back(vertex);
            normals.push_back(normal);
        }
    }

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < m; i++)
        {
            int idx1 = startIdx + j * (m + 1) + i;
            int idx2 = idx1 + 1;
            int idx3 = idx1 + (m + 1);
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

void generateGrid(int m, int n, const std::string &filename)
{
    generateCustomGrid(m, n, 0.0f, 2,
                       filename);
}

void generateCube(int m, int n, int o, const std::string &filename)
{
    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    std::vector<std::vector<int>> faces;

    int total = (m + 1) * (n + 1) * (o + 1);

    // Map to store unique surface vertices
    IntList verticeMap(total, -1);

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

    auto addFace = [&](int i1, int j1, int k1, int i2, int j2, int k2, int i3, int j3, int k3, int i4, int j4, int k4)
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

    // std::cerr << "Vertices"<<filename
    printerr("Actual Vertices : ", vertices.size(), " | Normals : ", normals.size(), " | Faces : ", faces.size());
    int exp = 8 + 4 * (m - 1 + n - 1 + o - 1) + 2 * ((m - 1) * (n - 1) + (n - 1) * (o - 1) + (o - 1) * (m - 1));
    printerr("Expect Vertices : ", exp, " | Normals : ", exp, " | Faces : ", 2 * (m * n + n * o + o * m));
    saveOBJ(filename, vertices, normals, faces);
}

void generateSphere(int m, int n, const std::string &filename, int axis, int direction)

{
    // axis : 0 for x and 2 for z
    assert(axis <= 2);
    assert(axis >= 0);
    std::vector<Vec3> vertices;
    std::vector<Vec3> normals;
    std::vector<std::vector<int>> faces;

    auto addFace = [&](std::initializer_list<int> indices)
    {
        std::vector<int> face(indices);
        std::reverse(face.begin(), face.end()); // Reverse before inserting
        faces.push_back(face);
    };

    Vec3 vec = {0.0f, 0.0f, 0.0f};
    // vec[axis] = 1;
    vec[axis] = direction;
    vertices.push_back(vec); // Top pole
    normals.push_back(vec);

    for (int i = 1; i < n; i++)
    {
        float phi = M_PI * i / n;
        for (int j = 0; j < m; j++)
        {
            float theta = 2 * M_PI * j / m;
            float x = sin(phi) * cos(theta);
            float y = sin(phi) * sin(theta);
            float z = cos(phi);
            // vertices.push_back({x, y, z});
            // normals.push_back({x, y, z});
            Vec3 v;
            if (axis == 0)
                v = {z, y, x}; // Swap X and Z
            else if (axis == 1)
                v = {x, z, y}; // Swap Y and Z
            else
                v = {x, y, z}; // Default (Z-axis)
            v[axis] = v[axis] * direction;

            vertices.push_back(v);
            normals.push_back(glm::normalize(v));
        }
    }

    vec[axis] = -vec[axis];
    // vec[axis] = vec[axis] * direction;

    vertices.push_back(vec); // Bottom pole
    normals.push_back(vec);

    for (int j = 0; j < m; j++)
    {
        addFace({(j + 1) % m + 1, j + 1, 0});
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
            addFace({idx1, idx2, idx4, idx3});
        }
    }

    int last = vertices.size() - 1;
    for (int j = 0; j < m; j++)
    {
        addFace({last, last - m + j, last - m + (j + 1) % m});
    }

    saveOBJ(filename, vertices, normals, faces);
}

// part 3

void MeshHalfEdge::loadObjfile(const std::string &filename, Vec3List &vertices,
                               Vec3List &normals,
                               FaceList &faces)
{

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
            // std::cout << "\tOne row of  vertices : " << v.x << " " << v.y << " " << v.z << std::endl;
            vertices.push_back(v);
        }
        else if (prefix == "vt")
        { // Texture coordinates
            continue;
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

                // Ignore texture information
                // if (!vtIdx.empty()) face.texIndices.push_back(std::stoi(vtIdx) - 1);

                // if (!vnIdx.empty())
                //     face.normalIndices.push_back(std::stoi(vnIdx) - 1);
            }
            faces.push_back(face);
        }
    }

    file.close();
    this->vertexPos = std::move(vertices);
    if (normals.size() < vertices.size())
    {
        normals.resize(vertices.size(), Vec3(0.0f, 0.0f, 0.0f)); // Extend with zero normals
    }
    this->vertexNormal = std::move(normals);
}

// part4
Vec3 MeshHalfEdge::computeFaceNormal(int x, int y, int z)
{

    Vec3 edge1 = vertexPos[z] - vertexPos[y];
    Vec3 edge2 = vertexPos[y] - vertexPos[x];

    Vec3 faceNormal = -glm::normalize(glm::cross(edge1, edge2));
    return faceNormal;
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
            u = halfEdge[h][HEAD];
            h = halfEdge[h][PAIR];

            assert(halfEdge[h][HEAD] == v);

            if (halfEdge[h][NEXT] == -1)
            {
                // Must be a boundary edge, so LEFT and NEXT both must be 0
                assert(halfEdge[h][LEFT] == -1);
                // printerr("This is a boundary edge:", h);
                break;
            }

            h = halfEdge[h][NEXT];

            w = halfEdge[h][HEAD];

            normal = normal + computeFaceNormal(u, v, w);
            // print(u, "->", v, "->", w);
        } while (h != start);
        vertexNormal[v] = glm::normalize(normal);
        // vertexNormal[v] = normal / (float)count;
    }
}

// part 5
void MeshHalfEdge::smoothen_step(float λ, int maxVertex)
{
    int start;
    int h;

    int v1;
    float count;
    if (maxVertex == -1)
        maxVertex = vertexPos.size();
    for (int v = 0; v < maxVertex; v++)
    {

        if (false)
            print("Smoothen vertex", v);
        Vec3 sum = Vec3(0.0f);

        count = 0;
        start = vertexHalfEdge[v]; // outgoing half edge
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

        assert(count > 0);

        // if (count == 0.0f)
        //     continue;
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

void MeshHalfEdge::smoothen_taubin(float λ, float u, int iterations, int maxVertex)
{
    // taubin
    for (int iter = 0; iter < iterations; iter++)
    {
        smoothen_step(λ, maxVertex);
        smoothen_step(u, maxVertex);
    }
}

void MeshHalfEdge::smoothen_laplacian(float λ, int iterations, int maxVertex)
{

    smoothen_taubin(λ, 0.0f, iterations, maxVertex);
}

// Function to generate a random displacement vector (Gaussian or Uniform)
Vec3 getRandomDisplacement(const std::string &noiseType, float param)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());

    if (noiseType == "gaussian")
    {
        std::normal_distribution<float> dist(0.0f, param); // Mean = 0, Std dev = sigma
        return Vec3(dist(gen), dist(gen), dist(gen));
    }
    else if (noiseType == "uniform")
    {
        std::uniform_real_distribution<float> dist(-param, param); // Range [-param, param]
        return Vec3(dist(gen), dist(gen), dist(gen));
    }
    else
    {
        throw std::invalid_argument("Invalid noise type! Use 'gaussian' or 'uniform'.");
    }
}

// Apply noise to all vertices in the mesh
void MeshHalfEdge::addNoise(const std::string &noiseType, float param)
{
    for (Vec3 &v : vertexPos)
    {
        v += getRandomDisplacement(noiseType, param);
    }
}

// part 6
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

// Part 8
void MeshHalfEdge::extrudeVertex(int v, float factor)
{
    if (v == -1)
        v = vertexPos.size() - 1;
    vertexPos[v] = vertexPos[v] + vertexNormal[v] * factor;
}

void MeshHalfEdge::extrudeNeighbors(int v, float factor)
{

    if (v == -1)
        v = vertexPos.size() - 1;

    int h, start, u;
    h = vertexHalfEdge[v];
    start = h;
    Vec3 normal = vertexNormal[v] * factor;

    do
    {
        u = halfEdge[h][HEAD];
        // use normal of v to extrude its neoughbour u
        vertexPos[u] = vertexPos[u] + normal;
        h = halfEdge[h][PAIR];
        h = halfEdge[h][NEXT];
    } while (h != start);
}

void MeshHalfEdge::extrudeCopyNeighbors(int v, float factor)
{

    if (v == -1)
        v = vertexPos.size() - 1;

    int h, start, u;
    h = vertexHalfEdge[v];
    start = h;
    Vec3 normal = vertexNormal[v] * factor;

    int x, y, z;
    IntList ring;
    do
    {
        ring.push_back(halfEdge[h][NEXT]);

        h = halfEdge[h][NEXT];
        h = halfEdge[h][NEXT];
        h = halfEdge[h][PAIR];
    } while (h != start);

    int next, pair, face;
    int first_next, V_old = vertexPos.size(), H_old = halfEdge.size();
    vertexHalfEdge[v] = H_old + 3;

    for (int i = 0; i < ring.size(); i++)
    {
        h = ring[i];
        u = halfEdge[h][HEAD];

        int V = vertexPos.size();
        int H = halfEdge.size();
        int F = FaceHalfEdge.size();
        vertexPos.push_back(vertexPos[u] + normal);
        vertexNormal.push_back(normal);

        FaceHalfEdge.push_back(H + 1);
        vertexHalfEdge.push_back(H + 1);

        next = halfEdge[h][NEXT];
        pair = halfEdge[next][PAIR];
        face = halfEdge[pair][LEFT];

        halfEdge[next][HEAD] = V;
        halfEdge[next][NEXT] = H - 4;

        // PAIR = 0,
        // NEXT,
        // HEAD,
        // LEFT

        halfEdge.push_back(std::array<int, 4>{H + 1, pair, V, face});
        halfEdge.push_back(std::array<int, 4>{H, H + 2, V + 1, F});
        halfEdge.push_back(std::array<int, 4>{H + 5, H + 3, v, F});
        halfEdge.push_back(std::array<int, 4>{H - 2, H + 1, V, F});

        if (i == 0)
        {
            // first
            first_next = next;
            // halfEdge[next][NEXT] = H - 4;
        }
        if (i == ring.size() - 1)
        {
            // last one
            halfEdge[first_next][NEXT] = H;
            halfEdge[H + 1][HEAD] = V_old;
            halfEdge[H + 2][PAIR] = H_old + 3;
        }
    }
}

void MeshHalfEdge::planarizeNeighbors(int v)
{
    // pull neighbours up along normal of v , aligning them along normal

    if (v == -1)
        v = vertexPos.size() - 1;

    int h, start, u;

    h = vertexHalfEdge[v];
    start = h;
    Vec3 normal = vertexNormal[v];

    do
    {
        u = halfEdge[h][HEAD];
        // use normal of v to extrude its neoughbour u

        float factor = glm::dot(vertexPos[v] - vertexPos[u], normal);
        vertexPos[u] = vertexPos[u] + normal * factor;
        h = halfEdge[h][PAIR];
        h = halfEdge[h][NEXT];
    } while (h != start);
}

void MeshHalfEdge::flatten(int v)
{

    if (v == -1)
        v = vertexPos.size() - 1;

    int count = 0;
    Vec3 pos = {0.0f, 0.0f, 0.0f};

    int h, start, u;

    h = vertexHalfEdge[v];
    start = h;

    do
    {

        printerr("Currently h = ", h);
        u = halfEdge[h][HEAD];

        // use normal of v to extrude its neoughbour u
        pos = pos + vertexPos[u];
        count++;
        h = halfEdge[h][PAIR];
        h = halfEdge[h][NEXT];
    } while (h != start);

    pos = pos / (float)count;
    vertexPos[v] = pos;
}
