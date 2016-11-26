#include "Sketch.h"
#include "Vertex.h"

void Vertex::register_parameters(Sketch *s) {
    s->add_parameter(X);
    s->add_parameter(Y);
};

std::pair<double, double> Vertex::rotate(std::shared_ptr<Vertex> const &origin, double const angle) const {
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