#include "Sketch.hpp"

bool Curve::on_segment(std::shared_ptr<Vertex> const &v) const {
    return on_segment(v->x(), v->y());
}

bool Curve::on_segment(std::shared_ptr<Vertex> const &v, std::shared_ptr<Vertex> const &origin, double const angle) const {
    double x, y;
    std::tie(x, y) = v->rotate(origin, angle);

    return on_segment(x, y);
}

bool Curve::on_manifold(std::shared_ptr<Vertex> const &v) const {
    return on_manifold(v->x(), v->y());
}

bool Curve::on_manifold(std::shared_ptr<Vertex> const &v, std::shared_ptr<Vertex> const &origin, double const angle) const {
    double x, y;
    std::tie(x, y) = v->rotate(origin, angle);

    return on_manifold(x, y);
}