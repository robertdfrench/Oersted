#include "Sketch.h"
#include "Vertex.h"

void Vertex::register_parameters(Sketch *s) const {
    s->add_parameter(std::const_pointer_cast<Variable>(X));
    s->add_parameter(std::const_pointer_cast<Variable>(Y));
};

double2 Vertex::rotate(std::shared_ptr<Vertex const> const &origin, double angle) const {
    double x = origin->x();
    double y = origin->y();

    double dx = this->x() - x;
    double dy = this->y() - y;

    double r = sqrt(dx * dx + dy * dy);
    double a = atan2(dy, dx) + angle * M_PI / 180.0;

    x += r * cos(a);
    y += r * sin(a);

    return double2{x,y};
}