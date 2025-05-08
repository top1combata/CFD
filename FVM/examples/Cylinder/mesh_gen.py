import sys
import os
from pathlib import Path

from boundary_condition import Boundaries, BoundaryType
from mesh_from_gmsh import mesh_from_gmsh
from mesh_to_cpp_string import mesh_to_cpp_string

from boundary_condition import wall_boundaries, outlet_boundaries, inlet_boundaries


INLET_VELOCITY = 0.0075

if __name__ == "__main__":

    directory = Path(__file__).resolve().parent
    with open(f"{directory}/cylinder.geo", "r") as file:
        gmsh_program = file.read()

    boundaries_map = {
        "wall": wall_boundaries(),
        "outlet": outlet_boundaries(0),
        "inlet": inlet_boundaries((INLET_VELOCITY, 0))
    }

    dir = os.path.dirname(os.path.abspath(__file__))
    with open(f"{dir}/cylinder_mesh.h", "w") as file:
        file.write(
            mesh_to_cpp_string(
                mesh_from_gmsh(gmsh_program, boundaries_map)
            )
        )
