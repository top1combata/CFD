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

    bool useAbsoluteVectorLength = true;
    auto velocityField = getVelocityField(m_solver, 0, useAbsoluteVectorLength);

    auto applyScale = [&velocityField](Scalar scale)
    {
        for (auto& arrow : velocityField)
        {
            arrow.setScale(scale);
        }
    };

    GuiWindow window(800, 600, "Post Processing");
    window.setView(getInitializeView(m_mesh));

    Index timePointIdx = 0;
    float scale = 1;

    while (window.isOpen())
    {
        window.pollEvent();
        window.clear({230, 230, 230});

        ImGui::Begin("Controls");

        if (m_solver.isTransient() && ImGui::SliderInt("Time Point", &timePointIdx, 0, m_solver.getTimePointAmount() - 1))
        {
            velocityField = getVelocityField(m_solver, timePointIdx, useAbsoluteVectorLength);
            applyScale(scale);
        }

        if (ImGui::Checkbox("Absolute vector size", &useAbsoluteVectorLength))
        {
            velocityField = getVelocityField(m_solver, timePointIdx, useAbsoluteVectorLength);
            scale = 1;
        }

        if (ImGui::SliderFloat("Vector Scale", &scale, 0, 10))
        {
            applyScale(scale);
        }

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
