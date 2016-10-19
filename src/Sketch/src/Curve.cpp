#include "Sketch.hpp"

Curve *Curve::split(Vertex *vnew, double s) {
    Curve *cnew = clone();
    *(vnew) = point(s);

    End = vnew;
    cnew->Start = vnew;

    return cnew;
}

bool Curve::on_segment(const Vertex *v) const {
    return on_segment(v->x(), v->y());
}

bool Curve::on_segment(const Vertex *v, const Vertex *origin, const double angle) const {
    double x, y;
    std::tie(x, y) = v->rotate(origin, angle);

    return on_segment(x, y);
}

bool Curve::on_manifold(const Vertex *v) const {
    return on_manifold(v->x(), v->y());
}

bool Curve::on_manifold(const Vertex *v, const Vertex *origin, const double angle) const {
    double x, y;
    std::tie(x, y) = v->rotate(origin, angle);

    return on_manifold(x, y);
}