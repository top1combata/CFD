#pragma once

#include "Utils/Types.h"
#include "ViewManager.h"

#include <SFML/Graphics.hpp>


class GuiWindow
{
public:

    GuiWindow(unsigned width, unsigned height, std::string const& title);
    ~GuiWindow();

    bool isOpen() const;
    void clear(sf::Color color = sf::Color::White);
    void renderGuiAndDisplay();
    void pollEvent();
    void render(sf::Drawable const& drawable);

    void setView(sf::View const& view);
    
private:

    sf::RenderWindow m_window;
    sf::Clock m_deltaClock;
    ViewManager m_viewManager;
};
