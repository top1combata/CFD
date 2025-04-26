#include "ViewManager.h"

#include <imgui-SFML.h>
#include <imgui.h>


ViewManager::ViewManager(sf::RenderWindow& window)
    : m_window(window)
    , m_prevWindowSize(window.getSize())
{}


void ViewManager::processEvent(sf::Event const& event)
{
    if (event.is<sf::Event::Resized>())
    {
        sf::Vector2u newSize = event.getIf<sf::Event::Resized>()->size;
        sf::View view = m_window.getView();
        view.setSize
        ({
            newSize.x / (float)m_prevWindowSize.x * view.getSize().x,
            newSize.y / (float)m_prevWindowSize.y * view.getSize().y
        });
        m_prevWindowSize = newSize;
        m_window.setView(view);

        return;
    }

    if (event.is<sf::Event::MouseWheelScrolled>())
    {
        auto const& scroll = event.getIf<sf::Event::MouseWheelScrolled>();
        if (scroll->wheel != sf::Mouse::Wheel::Vertical)
        {
            return;
        }

        float zoomFactor = 1 - zoomSensitivity * scroll->delta;
        sf::View view = m_window.getView();
        view.setSize(view.getSize() * zoomFactor);
        m_window.setView(view);

        return;
    }

    if (event.is<sf::Event::MouseButtonPressed>())
    {
        auto const& mouseEvent = event.getIf<sf::Event::MouseButtonPressed>();
        if (mouseEvent->button != sf::Mouse::Button::Left)
        {
            return;
        }

        m_isDragging = true;
        m_prevMousePosition = m_window.mapPixelToCoords(mouseEvent->position, m_window.getView());

        return;
    }

    if (event.is<sf::Event::MouseButtonReleased>())
    {
        auto const& mouseEvent = event.getIf<sf::Event::MouseButtonReleased>();
        if (mouseEvent->button != sf::Mouse::Button::Left)
        {
            return;
        }

        m_isDragging = false;

        return;
    }

    if (event.is<sf::Event::MouseMoved>())
    {
        // Prevent dragging while focusing on widgets
        if (!m_isDragging || ImGui::GetIO().WantCaptureMouse)
        {
            return;
        }

        auto const& mouseEvent = event.getIf<sf::Event::MouseMoved>();
        sf::Vector2f currMousePosition = m_window.mapPixelToCoords(mouseEvent->position, m_window.getView());
        sf::Vector2f delta = currMousePosition - m_prevMousePosition;

        sf::View view = m_window.getView();
        view.setCenter(view.getCenter() - delta);
        m_window.setView(view);
        m_prevMousePosition = m_window.mapPixelToCoords(mouseEvent->position, view);

        return;
    }
}
