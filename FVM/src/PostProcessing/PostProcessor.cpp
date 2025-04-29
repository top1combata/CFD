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
    GuiWindow window(800, 600, "Post Processing");
    window.setView(getInitializeView(m_mesh));

    auto meshCells = getMeshCells(m_mesh);
    auto meshFaces = getMeshFaces(m_mesh);
    List<Arrow> velocityField;

    enum PlotType
    {
        VELOCITY,
        VELOCITY_MAGNITUDE,
        PRESSURE
    };

    bool useAbsoluteVectorLength = true;
    Index plotType = VELOCITY;
    Index timePointIdx = 0;
    float scale = 1;

    auto updateEnitites = [&]()
    {
        auto applyScale = [&]()
        {
            for (auto& arrow : velocityField)
            {
                arrow.setScale(scale);
            }
        };

        switch (plotType)
        {
        case VELOCITY:
            velocityField = getVelocityField(m_solver, timePointIdx, useAbsoluteVectorLength);
            applyScale();
            break;

        case VELOCITY_MAGNITUDE:
            heatMapColor(meshCells, [this, &timePointIdx](Index cellIdx)
            {
                return m_solver.getVelocity(timePointIdx)(cellIdx).norm();
            });      
            break;

        case PRESSURE:
            heatMapColor(meshCells, [this, &timePointIdx](Index cellIdx)
            {
                return m_solver.getPressure(timePointIdx)(cellIdx);
            });
            break;
        }
    };

    auto renderEntities = [&]()
    {
        for (auto const& face : meshFaces)
        {
            window.render(face);
        }

        switch (plotType)
        {
        case VELOCITY:
            for (auto const& arrow : velocityField)
            {
                window.render(arrow);
            }
            break;

        case VELOCITY_MAGNITUDE:
        case PRESSURE:
            for (auto const& cell : meshCells)
            {
                window.render(cell);
            }
            break;
        }
    };

    updateEnitites();
    while (window.isOpen())
    {
        window.pollEvent();
        window.clear({230, 230, 230});

        ImGui::Begin("Controls");

        if (ImGui::Combo("Plot type", &plotType, "Velocity\0Velocity magnitude\0Pressure\0"))
        {
            updateEnitites();
        }

        if (m_solver.isTransient() && ImGui::SliderInt("Time point", &timePointIdx, 0, m_solver.getTimePointAmount() - 1))
        {
            updateEnitites();
        }

        if (plotType == VELOCITY)
        {
            if (ImGui::Checkbox("Absolute vector size", &useAbsoluteVectorLength))
            {
                updateEnitites();
            }

            if (ImGui::SliderFloat("Vector scale", &scale, 0, 10))
            {
                updateEnitites();
            }
        }

        ImGui::End();

        renderEntities();

        window.renderGuiAndDisplay();
    }
}
