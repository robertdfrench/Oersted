#include "Sketch.hpp"

std::pair<double, double> Vertex::rotate(const Vertex *origin, const double angle) const {
    double x = origin->x();
    double y = origin->y();

    double dx = this->x() - x;
    double dy = this->y() - y;

    double r = sqrt(dx * dx + dy * dy);
    double a = atan2(dy, dx) + angle * M_PI / 180.0;

    x += r * cos(a);
    y += r * sin(a);

    return std::make_pair(x, y);
}