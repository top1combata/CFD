import sys
from pathlib import Path

from boundary_condition import Boundaries, BoundaryType
from mesh_from_gmsh import mesh_from_gmsh

from boundary_condition import wall_boundaries, outlet_boundaries, inlet_boundaries


if __name__ == "__main__":

    directory = Path(__file__).resolve().parent
    with open(f"{directory}/cylinder.geo", "r") as file:
        gmsh_program = file.read()

    inlet_velocity = (0.15 if len(sys.argv) < 2  else float(sys.argv[1]))

    boundaries_map = {
        "wall": wall_boundaries(),
        "outlet": outlet_boundaries(0),
        "inlet": inlet_boundaries((inlet_velocity, 0))
    }

    with open("cylinder.msh", "w") as file:
        file.write(mesh_from_gmsh(gmsh_program, boundaries_map))
