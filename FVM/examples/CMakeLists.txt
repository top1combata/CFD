find_program(GMSH gmsh REQUIRED)
find_package(Python REQUIRED)

add_subdirectory(Cylinder)
add_subdirectory(HagenPoiseuille)
add_subdirectory(LidDrivenCavity)
add_subdirectory(Elbow)
