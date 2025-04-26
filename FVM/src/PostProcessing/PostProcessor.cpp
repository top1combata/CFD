#include "PostProcessor.h"
#include "Window/Window.h"
#include "Drawing/Drawing.h"

#include <imgui-SFML.h>
#include <imgui.h>


PostProcessor::PostProcessor(SolverBase const& solver)
    : m_mesh(solver.getMesh()), m_solver(solver)
{
    if (!m_mesh.is2D())
    {
        throw std::runtime_error("PostProcessor is only for 2D meshes !\n");
    }
}


void PostProcessor::show()
{
    auto meshCells = getMeshCells(m_mesh);
    auto meshFaces = getMeshFaces(m_mesh);
    auto velocityField = getVelocityField(m_solver);

    GuiWindow window(800, 600, "Post Processing");
    window.setView(getInitializeView(m_mesh));

    while (window.isOpen())
    {
        window.pollEvent();
        window.clear({230, 230, 230});

        ImGui::Begin("Controls");
        ImGui::End();

        for (auto const& face : meshFaces)
        {
            window.render(face);
        }
        for (auto const& arrow : velocityField)
        {
            window.render(arrow);
        }

        window.renderGuiAndDisplay();
    }
}
