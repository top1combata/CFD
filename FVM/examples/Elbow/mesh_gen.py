import sys
from pathlib import Path

from boundary_condition import Boundaries, BoundaryType
from mesh_from_gmsh import mesh_from_gmsh

from boundary_condition import wall_boundaries, outlet_boundaries, inlet_boundaries


if __name__ == "__main__":

    directory = Path(__file__).resolve().parent
    with open(f"{directory}/elbow.geo", "r") as file:
        gmsh_program = file.read()

    main_inlet_velocity = (0.1 if len(sys.argv) < 2  else float(sys.argv[1]))
    secondary_inlet_velocity = (0.5 if len(sys.argv) < 3 else float(sys.argv[2]))

    boundaries_map = {
        "wall": wall_boundaries(),
        "outlet": outlet_boundaries(0),
        "main-inlet": inlet_boundaries((main_inlet_velocity, 0)),
        "secondary-inlet": inlet_boundaries((0, secondary_inlet_velocity))
    }

    with open("elbow.msh", "w") as file:
        file.write(mesh_from_gmsh(gmsh_program, boundaries_map))
