#include "MeshDrawer.h"
#include "Window/Window.h"
#include "Drawing/Drawing.h"


MeshDrawer::MeshDrawer(MeshBase const& mesh)
    : m_mesh(mesh)
{
    if (!m_mesh.is2D())
    {
        throw std::runtime_error("PostProcessor is only for 2D meshes !\n");
    }
}


void MeshDrawer::show()
{
    GuiWindow window(800, 600, "Mesh Drawer");
    window.setView(getInitializeView(m_mesh));

    auto meshCells = getMeshCells(m_mesh);
    auto meshFaces = getMeshFaces(m_mesh);

    while (window.isOpen())
    {
        window.pollEvent();
        window.clear({230, 230, 230});

        for (auto const& cell : meshCells)
        {
            window.render(cell);
        }

        for (auto const& face : meshFaces)
        {
            window.render(face);
        }

        window.renderGuiAndDisplay();
    }
}
