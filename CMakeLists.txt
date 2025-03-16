cmake_minimum_required(VERSION 3.0)

project(a2)
# set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD 20) #Needed to used span
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(glm REQUIRED)
find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)

add_library(viewer src/hw.cpp src/viewer.cpp deps/src/gl.c)
target_include_directories(viewer PUBLIC deps/include)
target_link_libraries(viewer glm::glm OpenGL::GL SDL2::SDL2)

add_executable(example src/example.cpp)
target_link_libraries(example viewer)


add_executable(e3 src/e3.cpp src/mesh.cpp)
target_link_libraries(e3 viewer)

add_executable(e5 src/e5.cpp src/mesh.cpp)
target_link_libraries(e5 viewer)

add_executable(e6 src/e6.cpp src/mesh.cpp)
target_link_libraries(e6 viewer)

add_executable(e7 src/e7.cpp src/mesh.cpp)
target_link_libraries(e7 viewer)

add_executable(e8 src/e8.cpp src/mesh.cpp)
target_link_libraries(e8 viewer)

add_executable(e9 src/e9.cpp src/mesh.cpp)
target_link_libraries(e9 viewer)
