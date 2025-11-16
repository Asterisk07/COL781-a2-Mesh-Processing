# COL781 Assignment 2: Mesh Processing

### This top I modelled , along with the reference image:

<img width="1920" height="1020" alt="7_4" src="https://github.com/user-attachments/assets/4c06e84e-811a-4590-9834-149481e3a7a9" />
<img width="655" height="511" alt="7_3" src="https://github.com/user-attachments/assets/f83700a4-d6cb-4c6c-89a0-383b08664b44" />
<img width="644" height="505" alt="7_6" src="https://github.com/user-attachments/assets/fc00ce9a-8bb7-417c-b713-d00ac332dd0d" />

### Normal recomputation for the given meshes:

<img width="643" height="510" alt="1_spot (1)" src="https://github.com/user-attachments/assets/58a42755-f72e-4cf5-addf-2964f82afad5" />
<img width="643" height="506" alt="1_bunny (1)" src="https://github.com/user-attachments/assets/6625bb86-e3d9-477a-9e7a-ba22f5c07565" />


Make sure that [glm](https://github.com/g-truc/glm) and [SDL2](https://www.libsdl.org/) are installed. Ideally, these should be installed by your package manager rather than manually (at least, if you are on Linux or Mac).

Then compile the code using the standard CMake procedure:

- The first time, run `cmake -B build` from the project root to create a `build/` directory and initialize a build system there.
- Then, every time you want to compile the code, run `cmake --build build` (again from the project root). Then the example programs will be created under `build/`.
