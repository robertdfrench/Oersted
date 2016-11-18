#include "Sketch.hpp"

std::pair<double, double> Vertex::rotate(std::shared_ptr<Vertex> const &origin, double const angle) const {
    double x = origin->x();
    double y = origin->y();

    double dx = this->x() - x;
    double dy = this->y() - y;

    double r = std::sqrt(dx * dx + dy * dy);
    double a = std::atan2(dy, dx) + angle * M_PI / 180.0;

    x += r * std::cos(a);
    y += r * std::sin(a);

    return std::make_pair(x, y);
}