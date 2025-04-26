#include "Window.h"

#include <imgui-SFML.h>
#include <imgui.h>


sf::ContextSettings getContextSettings()
{
    sf::ContextSettings settings;
    settings.antiAliasingLevel = 2;
    return settings;
}

GuiWindow::GuiWindow(unsigned width, unsigned height, std::string const& title)
    : m_window(sf::VideoMode({width, height}), title, sf::State::Windowed, getContextSettings())
    , m_viewManager(m_window)
{
    if (!ImGui::SFML::Init(m_window, true))
    {
        throw std::runtime_error("failed to initialize gui");
    }
    m_window.setFramerateLimit(60);
    m_window.setVerticalSyncEnabled(true);
}


GuiWindow::~GuiWindow()
{
    ImGui::SFML::Shutdown();
}


bool GuiWindow::isOpen() const
{
    return m_window.isOpen();
}


void GuiWindow::clear(sf::Color color)
{
    m_window.clear(color);
}


void GuiWindow::renderGuiAndDisplay()
{
    ImGui::SFML::Render(m_window);
    m_window.display();
}


void GuiWindow::render(sf::Drawable const& drawable)
{
    m_window.draw(drawable);
}


void GuiWindow::pollEvent()
{
    while (const std::optional event = m_window.pollEvent())
    {
        if (event->is<sf::Event::Closed>())
        {
            m_window.close();
        }

        ImGui::SFML::ProcessEvent(m_window, *event);
        m_viewManager.processEvent(*event);
    }
    // Some ImGui magic
    ImGui::SFML::Update(m_window, m_deltaClock.restart());
}


void GuiWindow::setView(sf::View const& view)
{
    m_window.setView(view);
}
