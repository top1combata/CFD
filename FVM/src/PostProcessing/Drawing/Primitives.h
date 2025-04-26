#pragma once

#include "Utils/Types.h"

#include <SFML/Graphics.hpp>


class Line : public sf::Drawable
{
public:

    Line() = default;
    Line(Vector, Vector);

    void setColor(sf::Color color);

    void draw(sf::RenderTarget& target, sf::RenderStates states) const override;

private:

    sf::Vertex m_vertices[2];
};


class Arrow : public sf::Drawable
{
public:

    Arrow(Vector center, Vector direction);

    // Scale is from 0 ot 1
    void setScale(Scalar scale);
    void setDirection(Vector direction);

    void draw(sf::RenderTarget& target, sf::RenderStates states) const override;

private:

    static constexpr Scalar HEAD_ANGLE = 0.25;
    static constexpr Scalar RELATIVE_HEAD_SIZE = 0.2;

    Vector m_center;
    Vector m_direction;
    Scalar m_scale = 1;
    Line m_lines[3];

    void updateLines();
};
