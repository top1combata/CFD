#pragma once

#pragma once

#include <SFML/Graphics/View.hpp>
#include <SFML/Graphics/RenderWindow.hpp>


class ViewManager
{
public:

    ViewManager(sf::RenderWindow& window);

    void processEvent(sf::Event const& event);

private:

    static constexpr float zoomSensitivity = 0.03f;

    sf::RenderWindow& m_window;
    bool m_isDragging = false;
    sf::Vector2f m_prevMousePosition;
    sf::Vector2u m_prevWindowSize;
};
