def mesh_to_cpp_string(mesh_content: str) -> str:
    
    res = ""
    for line in mesh_content.split("\n"):
        res += "\"" + line + "\\n\"\n"
    
    return res
