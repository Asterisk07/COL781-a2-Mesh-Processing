#ifndef VIEWER_HPP
#define VIEWER_HPP

#include "hw.hpp"

namespace COL781 {
    namespace Viewer {

        class Camera {
        public:
            glm::vec3 position;
            glm::vec3 front;
            glm::vec3 up; 
            glm::vec3 lookAt;
            glm::mat4 viewMatrix;

            float cameraSpeed, yaw, pitch, lastX, lastY, fov, aspect;
            bool firstMouse;
            void initialize(float aspect);
            glm::mat4 getViewMatrix();
            glm::mat4 getProjectionMatrix();
            glm::vec3 getViewDir();
            glm::vec3 getRightVector();

            void setCameraView(glm::vec3 position_vector, glm::vec3 lookat_vector, glm::vec3 up_vector);
            void updateViewMatrix();
        };

        class Viewer {
        public:
            bool initialize(const std::string &title, int width, int height);
            void setMesh(int nv, int nt, int ne, const glm::vec3* vertices, const glm::ivec3* triangles, const glm::ivec2* edges, const glm::vec3* normals = nullptr);
            void view();
        private:
            COL781::OpenGL::Rasterizer r;
            COL781::OpenGL::ShaderProgram program;
            COL781::OpenGL::Object object;
            COL781::OpenGL::Object wireframe;
            glm::mat4 stagetransform;
            Camera camera;
        };

    }
}

#endif
