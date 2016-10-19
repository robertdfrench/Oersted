#include "Sketch.hpp"

Vertex LineSegment::point(double s) const {
    const double x0 = start()->x();
    const double y0 = start()->y();
    const double x1 = end()->x();
    const double y1 = end()->y();

    return Vertex{x0 * (1.0 - s) + x1 * s, y0 * (1.0 - s) + y1 * s};
}

Vertex LineSegment::tangent(double s, bool orientation) const {
    double dx = end()->x() - start()->x();
    double dy = end()->y() - start()->y();
    double dl = sqrt(dx * dx + dy * dy);

    dx /= dl;
    dy /= dl;

    if (orientation) {
        return Vertex{dx, dy};
    } else {
        return Vertex{-dx, -dy};
    }
}

double LineSegment::a(double s, bool orientation) const {
    if (orientation) {
        return atan2(end()->y() - start()->y(), end()->x() - start()->x());
    } else {
        return atan2(start()->y() - end()->y(), start()->x() - end()->x());
    }
};

double LineSegment::supremum() const {
    return std::fmax(start()->hypot(), end()->hypot());
};

bool LineSegment::on_manifold(const double x, const double y) const {
    double dxv = (x - 0.5 * (start()->x() + end()->x()));
    double dyv = (y - 0.5 * (start()->y() + end()->y()));
    double dv = sqrt(dxv * dxv + dyv * dyv);

    double dxl = start()->x() - end()->x();
    double dyl = start()->y() - end()->y();
    double dl = sqrt(dxl * dxl + dyl * dyl);

    double dot = abs(dxv * dxl + dyv * dyl);
    double norm = (dv * dl);
    double tol = fmax(dv, dl) * FLT_EPSILON;

    if (abs(dot - norm) < tol) {
        return true;
    } else {
        return false;
    }
}

bool LineSegment::on_segment(const double x, const double y) const {
    double dx0 = x - start()->x();
    double dy0 = y - start()->y();
    double dl0 = sqrt(dx0 * dx0 + dy0 * dy0);

    double dx1 = x - end()->x();
    double dy1 = y - end()->y();
    double dl1 = sqrt(dx1 * dx1 + dy1 * dy1);

    double dxl = start()->x() - end()->x();
    double dyl = start()->y() - end()->y();
    double dll = sqrt(dxl * dxl + dyl * dyl);

    double diff = abs(dl0 + dl1 - dll);
    double tol = dll * FLT_EPSILON;

    if (diff < tol) {
        // If the Vertex and two endpoints of the LineSegment form a degenerate triangle
        // Where the length of the legs attached to the Vertex equal the length of the LineSegment
        // Then the Vertex is on the LineSegment between the two endpoints
        return true;
    } else {
        return false;
    }
}

bool LineSegment::is_identical(const Curve *c) const {
    const LineSegment *l = dynamic_cast<const LineSegment *>(c);

    if (l == nullptr) {
        return false;
    } else {
        return is_identical(l->start()->x(), l->start()->y(), l->end()->x(), l->end()->y());
    }
}

bool LineSegment::is_identical(const Curve *c, const Vertex *origin, const double angle) const {
    const LineSegment *l = dynamic_cast<const LineSegment *>(c);

    if (l == nullptr) {
        return false;
    } else {
        double xs, ys, xe, ye;

        std::tie(xs, ys) = c->start()->rotate(origin, angle);
        std::tie(xe, ye) = c->end()->rotate(origin, angle);

        return is_identical(xs, ys, xe, ye);
    }
}

bool LineSegment::is_identical(const double x0, const double y0, const double x1, const double y1) const {
    double xs = start()->x();
    double ys = start()->y();

    double xe = end()->x();
    double ye = end()->y();

    double tol =
            FLT_EPSILON * std::fmax(abs(xs - xe), abs(ys - ye)); // #TODO: L1 norm is more efficient tolerance strategy

    return (abs(xs - x0) < tol && abs(ys - y0) < tol && abs(xe - x1) < tol && abs(ye - y1) < tol)
           || (abs(xs - x1) < tol && abs(ys - y1) < tol && abs(xe - x0) < tol && abs(ye - y0) < tol);
}

bool LineSegment::is_coincident(const Curve *c) const {
    const LineSegment *l = dynamic_cast<const LineSegment *>(c);

    if (l == nullptr) {
        return false;
    } else {
        if (on_manifold(c->start()) && on_manifold(c->end())) {
            return true;
        } else {
            return false;
        }
    }
}

// bool LineSegment::is_coincident(const Curve* c, const Vertex* origin, const double_t angle) {return on_manifold(c->start(), origin, angle) && on_manifold(c->end(), origin, angle)};

double LineSegment::length() const { return hypot(end()->x() - start()->x(), end()->y() - start()->y()); };

void LineSegment::replace_verticies(std::vector<Vertex *> oldv, std::vector<Vertex *> newv) {
    auto i = std::find(oldv.begin(), oldv.end(), Start);
    if (i != oldv.end()) {
        size_t j = i - oldv.begin();
        Start = newv[j];
    }

    i = std::find(oldv.begin(), oldv.end(), End);
    if (i != oldv.end()) {
        size_t j = i - oldv.begin();
        End = newv[j];
    }
}