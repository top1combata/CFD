from boundary_condition import Boundaries

from typing import Dict
from collections import defaultdict
from dataclasses import dataclass
import subprocess
import os


def mesh_from_gmsh(gmsh_program: str, boundaries_map: Dict[str, Boundaries]) -> str:
    
    SCRIPT_FILE_NAME = ".tmpscipt.geo"
    MESH_FILE_NAME = ".tmpmesh.msh"

    with open(SCRIPT_FILE_NAME, "w") as file:
        file.write(gmsh_program)

    # Generating mesh in msh2 format
    subprocess.run(["gmsh", "-2", SCRIPT_FILE_NAME, "-format",  "msh2", "-o", MESH_FILE_NAME])

    with open(MESH_FILE_NAME, "r") as file:
        lines = file.readlines()

    # Clear tmp files
    os.remove(SCRIPT_FILE_NAME)
    os.remove(MESH_FILE_NAME)

    return create_mesh(lines, boundaries_map)


@dataclass
class Node:
    x: float
    y: float
    z: float


def create_mesh(lines, boundaries_map: Dict[str, Boundaries]) -> str:
    lines_iter = iter(lines)
    line = next(lines_iter)

    phys_tag_to_name = {}
    # Physical names
    while "PhysicalNames" not in line:
        line = next(lines_iter)

    physical_names_count = int(next(lines_iter))
    for _ in range(physical_names_count):
        line = next(lines_iter)
        tokens = line.split(" ")
        phys_tag_to_name[int(tokens[1])] = tokens[2][1:-2]

    while "Nodes" not in line:
        line = next(lines_iter)

    nodes_count = int(next(lines_iter))
    nodes = []
    for _ in range(nodes_count):
        line = next(lines_iter)
        tokens = line.split(" ")
        nodes.append(Node(
            x=float(tokens[1]), 
            y=float(tokens[2]), 
            z=float(tokens[3])
        ))

    while "Elements" not in line:
        line = next(lines_iter)

    elements_count = int(next(lines_iter))
    physical_name_to_faces = defaultdict(list)
    faces = set()
    cells = []

    def sorted_pair(idx1, idx2):
        return (min(idx1, idx2), max(idx1, idx2))

    def add_face(node1, node2):
        nonlocal faces
        faces.add(sorted_pair(node1, node2))


    for _ in range(elements_count):
        line = next(lines_iter)
        tokens = line.split(" ")
        element_type = int(tokens[1])
        num_tags = int(tokens[2])
        physical_group = int(tokens[3])

        tokens = tokens[(3 + num_tags):]

        # Boundary Face
        if element_type == 1:
            node1 = int(tokens[0]) - 1
            node2 = int(tokens[1]) - 1
            add_face(node1, node2)
            physical_name_to_faces[phys_tag_to_name[physical_group]].append(sorted_pair(node1, node2))

        # Interior cell
        if element_type in [2, 3]:
            num_nodes = element_type + 1
            tokens.append(tokens[0])
            cell_nodes = []

            for i in range(num_nodes):
                node1 = int(tokens[i]) - 1
                node2 = int(tokens[i+1]) - 1

                add_face(node1, node2)
                cell_nodes.append(sorted_pair(node1, node2))

            cells.append(cell_nodes)
    
    face_amount = 0
    face_to_idx = {}
    for face in faces:
        face_to_idx[face] = face_amount
        face_amount += 1

    mesh_file_content = ""

    for node in nodes:
        mesh_file_content += f"v {node.x} {node.y}\n"

    for face in faces:
        mesh_file_content += f"f {face[0]} {face[1]}\n"

    for cell in cells:
        mesh_file_content += "c "
        for face in cell:
            mesh_file_content += f"{face_to_idx[face]} "
        mesh_file_content += "\n"

    def write_boundaries(face_idx, boundaries: Boundaries):
        nonlocal mesh_file_content

        p_boundary = boundaries.p_boundary
        u_boundary = boundaries.u_boundary

        mesh_file_content += f"b {face_idx} p {p_boundary.type.name} {p_boundary.value}\n"
        mesh_file_content += f"b {face_idx} U {u_boundary.type.name} {u_boundary.value[0]} {u_boundary.value[1]}\n"


    for physical_tag_name, faces in physical_name_to_faces.items():
        boundaries = boundaries_map[physical_tag_name]

        for face in faces:
            write_boundaries(face_to_idx[face], boundaries)

    return mesh_file_content
