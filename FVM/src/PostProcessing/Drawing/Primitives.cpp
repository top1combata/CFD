#include "Primitives.h"


Line::Line(Vector begin, Vector end)
{
    m_vertices[0].position = sf::Vector2f(begin.x(), begin.y());
    m_vertices[1].position = sf::Vector2f(end.x(), end.y());
    
    setColor(sf::Color::Black);
}


void Line::setColor(sf::Color color)
{
    m_vertices[0].color = color;
    m_vertices[1].color = color;
}


void Line::draw(sf::RenderTarget& target, sf::RenderStates states) const
{
    target.draw(m_vertices, 2, sf::PrimitiveType::Lines, states);
}


Arrow::Arrow(Vector center, Vector direction)
    : m_center(center), m_direction(direction)
{
    updateLines();
}


void Arrow::setScale(Scalar scale)
{
    m_scale = scale;
    updateLines();
}


void Arrow::setDirection(Vector direction)
{
    m_direction = direction;
    updateLines();
}


void Arrow::draw(sf::RenderTarget& target, sf::RenderStates states) const
{
    for (Line const& line : m_lines)
    {
        line.draw(target, states);
    }
}


void Arrow::updateLines()
{
    const Scalar angle = HEAD_ANGLE;

    Vector head = m_center + m_scale * m_direction;

    Vector delta =
    {
        cos(angle) * m_direction.x() - sin(angle) * m_direction.y(),
        sin(angle) * m_direction.x() + cos(angle) * m_direction.y(),
        0
    };
    delta *= m_scale * RELATIVE_HEAD_SIZE;
    Vector leftEdge = head - delta;

    delta = 
    {
        cos(angle) * m_direction.x() + sin(angle) * m_direction.y(),
        -sin(angle) * m_direction.x() + cos(angle) * m_direction.y(),
        0
    };
    delta *= m_scale * RELATIVE_HEAD_SIZE;
    Vector rightEdge = head - delta;

    m_lines[0] = Line(m_center, head);
    m_lines[1] = Line(head, leftEdge);
    m_lines[2] = Line(head, rightEdge);
}
