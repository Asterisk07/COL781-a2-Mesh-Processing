#include "viewer.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
namespace COL781 {
    namespace Viewer {

        namespace GL = COL781::OpenGL;

        void Camera::initialize(float aspect) {
            firstMouse = true;
            yaw   = -90.0f;    
            pitch =  0.0f;
            lastX =  800.0f / 2.0;
            lastY =  600.0 / 2.0;
            fov   =  60.0f;

            this->aspect = aspect;

            position = glm::vec3(0.0f, 0.0f,  1.5f);
            lookAt = glm::vec3(0.0f, 0.0f, 0.0f);
            up = glm::vec3(0.0f, 1.0f,  0.0f);

            updateViewMatrix();
        }

        glm::mat4 Camera::getViewMatrix() {
            return viewMatrix;
        }

        void Camera::updateViewMatrix() {
            viewMatrix = glm::lookAt(position, lookAt, up);
        }

        glm::mat4 Camera::getProjectionMatrix() {
            return glm::perspective(glm::radians(fov), aspect, 0.1f, 100.0f);
        }
            glm::vec3 getRightVector();

        glm::vec3 Camera:: getViewDir() {
            return -glm::transpose(viewMatrix)[2];
        }

        glm::vec3 Camera::getRightVector() {
            return glm::transpose(viewMatrix)[0];
        }

        void Camera::setCameraView(glm::vec3 position_vector, glm::vec3 lookat_vector, glm::vec3 up_vector) {
            position = std::move(position_vector);
            lookAt = std::move(lookat_vector);
            up = std::move(up_vector);

            viewMatrix = glm::lookAt(position, lookAt, up);
        }

        bool Viewer::initialize(const std::string &title, int width, int height) {
            if (!r.initialize(title.c_str(), width, height))
                return false;
            program = r.createShaderProgram(
                r.vsBlinnPhong(),
                r.fsBlinnPhong()
            );
            r.useShaderProgram(program);
            object = r.createObject();
            wireframe = r.createObject();
            r.enableDepthTest();
            camera.initialize((float)width/(float)height);
            return true;
        }

        glm::mat4 calculateModelMatrix(const glm::vec3* vertices, int nv) {
            if (nv <= 0) {
                return glm::mat4(1.0f);
            }
            glm::vec3 min = vertices[0];
            glm::vec3 max = vertices[0];
            for (int i = 1; i < nv; ++i) {
                min = glm::min(min, vertices[i]);
                max = glm::max(max, vertices[i]);
            }
            glm::vec3 center = (min + max) * 0.5f;
            glm::vec3 size = max - min;
            float max_dimension = glm::max(glm::max(size.x, size.y), size.z);
            if (max_dimension <= 0.0f) {
                max_dimension = 1.0f;
            }
            float scale_factor = 1.0f / max_dimension;
            glm::mat4 translation = glm::translate(glm::mat4(1.0f), -center);
            glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(scale_factor));
            return scaling * translation;
        }

        void Viewer::setMesh(int nv, int nt, int ne, const glm::vec3* vertices, const glm::ivec3* triangles, const glm::ivec2* edges, const glm::vec3* normals) {
            if(normals == nullptr) {
                glm::vec3* normalsz = new glm::vec3[nv];
                for(int i = 0; i < nv; i++) {
                    normalsz[i] = glm::vec3(0.0, 0.0, 0.0);
                }
                normals = normalsz;
            }
            r.setVertexAttribs(object, 0, nv, vertices);
            r.setVertexAttribs(object, 1, nv, normals);
            r.setTriangleIndices(object, nt, triangles);
            r.setVertexAttribs(wireframe, 0, nv, vertices);
            r.setVertexAttribs(wireframe, 1, nv, normals);
            r.setEdgeIndices(wireframe, ne, edges);
            stagetransform = calculateModelMatrix(vertices, nv);
        }

        void Viewer::view() {
            // The transformation matrix.
            glm::mat4 model = stagetransform;
            // glm::mat4 model = glm::mat4(1.0f);
            glm::mat4 view;    
            glm::mat4 projection = camera.getProjectionMatrix();

            float deltaAngleX = 2.0 * 3.14 / 800.0;
            float deltaAngleY = 3.14 / 600.0;

            int lastxPos, lastyPos, xPos, yPos;

            SDL_GetMouseState(&lastxPos, &lastyPos);

            while (!r.shouldQuit()) {
                r.clear(glm::vec4(1.0, 1.0, 1.0, 1.0));

                camera.updateViewMatrix();

                Uint32 buttonState = SDL_GetMouseState(&xPos, &yPos);
                if( buttonState & SDL_BUTTON(SDL_BUTTON_LEFT) ) {
                    glm::vec4 pivot = glm::vec4(camera.lookAt.x, camera.lookAt.y, camera.lookAt.z, 1.0f);
                    glm::vec4 position = glm::vec4(camera.position.x, camera.position.y, camera.position.z, 1.0f);

                    float xAngle = (float)(lastxPos - xPos) * deltaAngleX;
                    float yAngle = (float)(lastyPos - yPos) * deltaAngleY;

                    float cosAngle = dot(camera.getViewDir(), camera.up);

                    if(cosAngle * signbit(deltaAngleY) > 0.99f)
                        deltaAngleY = 0.0f;

                    glm::mat4 rotationMatX(1.0f);
                    rotationMatX = glm::rotate(rotationMatX, xAngle, camera.up);
                    position = (rotationMatX * (position - pivot)) + pivot;

                    glm::mat4 rotationMatY(1.0f);
                    rotationMatY = glm::rotate(rotationMatY, yAngle, camera.getRightVector());
                    glm::vec3 finalPosition = (rotationMatY * (position - pivot)) + pivot;
                    camera.position = finalPosition;
                    camera.updateViewMatrix();
                }

                buttonState = SDL_GetMouseState(&xPos, &yPos);
                if( buttonState & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
                    // Update camera parameters

                    float deltaY =  (float)(lastyPos - yPos) * 0.01f;
                    glm::mat4 dollyTransform = glm::mat4(1.0f);
                    dollyTransform = glm::translate(dollyTransform, normalize(camera.lookAt - camera.position) * deltaY);
                    glm::vec3 newCameraPosition = dollyTransform * glm::vec4(camera.position, 1.0f);
                    float newCameraFov = 2 * glm::atan(600.0f / (2 * deltaY)); // TODO Ask
                    
                    if(signbit(newCameraPosition.z) == signbit(camera.position.z)) {
                        camera.position = newCameraPosition;
                        camera.fov = newCameraFov; // TODO Ask
                        }
                }

                lastxPos = xPos;
                lastyPos = yPos;

                view = camera.getViewMatrix();
                
                r.setUniform(program, "model", model);
                r.setUniform(program, "view", view);
                r.setUniform(program, "projection", projection);
                r.setUniform(program, "lightPos", camera.position);
                r.setUniform(program, "viewPos", camera.position);
                r.setUniform(program, "lightColor", glm::vec3(1.0f, 1.0f, 1.0f));

                // r.setupFilledFaces();
                glm::vec3 red(1.0f, 0.0f, 0.0f);
                glm::vec3 green(0.0f, 1.0f, 0.0f);
                glm::vec3 blue(0.0f, 0.0f, 1.0f);
                glm::vec3 white(1.0f, 1.0f, 1.0f);
                r.setUniform(program, "ambientColor", 0.25f*white);
                r.setUniform(program, "intdiffuseColor", 0.75f*red);
                r.setUniform(program, "extdiffuseColor", 0.75f*blue);
                r.setUniform(program, "specularColor", 0.8f*white);
                r.setUniform(program, "phongExponent", 100.f);
                r.drawObject(object);

                r.setupWireFrame();
                glm::vec3 black(0.0f, 0.0f, 0.0f);
                r.setUniform(program, "ambientColor", black);
                r.setUniform(program, "intdiffuseColor", black);
                r.setUniform(program, "extdiffuseColor", black);
                r.setUniform(program, "specularColor", black);
                r.setUniform(program, "phongExponent", 0.f);
                r.drawEdges(wireframe);
                r.show();
            }
        }

    }
}
