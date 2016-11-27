#include "Curve.h"

bool Curve::on_segment(std::shared_ptr<Vertex const> const &v) const {
    return on_segment(v->x(), v->y());
}

bool Curve::on_segment(std::shared_ptr<Vertex const> const &v, std::shared_ptr<Vertex const> const &origin, double angle) const {
    double2 p = v->rotate(origin, angle);

    return on_segment(p.X, p.Y);
}

bool Curve::on_manifold(std::shared_ptr<Vertex const> const &v) const {
    return on_manifold(v->x(), v->y());
}

bool Curve::on_manifold(std::shared_ptr<Vertex const> const &v, std::shared_ptr<Vertex const> const &origin, double angle) const {
    double2 p = v->rotate(origin, angle);

    return on_manifold(p.X, p.Y);
}