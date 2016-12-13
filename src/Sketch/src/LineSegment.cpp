#include "LineSegment.h"
#include "doublen.h"

double2 LineSegment::point(double s) const {
    double x0 = start()->x();
    double y0 = start()->y();
    double x1 = end()->x();
    double y1 = end()->y();

    return double2{x0 * (1.0 - s) + x1 * s, y0 * (1.0 - s) + y1 * s};
}

double2 LineSegment::tangent(double s, bool orientation) const {
    double dx = end()->x() - start()->x();
    double dy = end()->y() - start()->y();
    double dl = sqrt(dx * dx + dy * dy);

    dx /= dl;
    dy /= dl;

    if (orientation) {
        return double2{dx, dy};
    } else {
        return double2{-dx, -dy};
    }
}

double LineSegment::a(double s, bool orientation) const {
    if (orientation) {
        return atan2(end()->y() - start()->y(), end()->x() - start()->x());
    } else {
        return atan2(start()->y() - end()->y(), start()->x() - end()->x());
    }
};

double2 LineSegment::supremum() const {
    double xs = start()->x();
    double ys = start()->y();
    double ls = sqrt(xs * xs + ys * ys);

    double xe = end()->x();
    double ye = end()->y();
    double le = sqrt(xe * xe + ye * ye);

    double dx = xe - xs;
    double dy = ye - ys;
    double dl = sqrt(dx * dx + dy * dy);
    dx /= dl;
    dy /= dl;

    double sup, cross;
    if (ls >= le) {
        sup = ls;
        cross = abs(xs * dy - ys * dx) / ls;
    } else {
        sup = le;
        cross = abs(xe * dy - ye * dx) / le;
    }

    return double2{sup, cross};
};

bool LineSegment::on_manifold(double x, double y) const {
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

bool LineSegment::on_segment(double x, double y) const {
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

MatchOrientation LineSegment::is_identical(std::shared_ptr<Curve const> const &c) const {
    std::shared_ptr<LineSegment const> l = std::dynamic_pointer_cast<LineSegment const>(c);

    if (l.get() == nullptr) {
        return MatchOrientation::None;
    } else {
        return is_identical(l->start()->x(), l->start()->y(), l->end()->x(), l->end()->y());
    }
}

MatchOrientation LineSegment::is_identical(std::shared_ptr<Curve const> const &c, std::shared_ptr<Vertex const> const &origin, double angle) const {
    std::shared_ptr<LineSegment const> l = std::dynamic_pointer_cast<LineSegment const>(c);

    if (l.get() == nullptr) {
        return MatchOrientation::None;
    } else {
        double2 ps = c->start()->rotate(origin, angle);
        double2 pe = c->end()->rotate(origin, angle);

        return is_identical(ps.X, ps.Y, pe.X, pe.Y);
    }
}

MatchOrientation LineSegment::is_identical(double x0, double y0, double x1, double y1) const {
    double xs = start()->x();
    double ys = start()->y();

    double xe = end()->x();
    double ye = end()->y();

    double tol = FLT_EPSILON * std::fmax(abs(xs - xe), abs(ys - ye));

    if (abs(xs - x0) < tol && abs(ys - y0) < tol && abs(xe - x1) < tol && abs(ye - y1) < tol) {
        return MatchOrientation::Forward;
    } else if (abs(xs - x1) < tol && abs(ys - y1) < tol && abs(xe - x0) < tol && abs(ye - y0) < tol) {
        return MatchOrientation::Reverse;
    } else {
        return MatchOrientation::None;
    }
}

bool LineSegment::is_coincident(std::shared_ptr<Curve const> const &c) const {
    std::shared_ptr<LineSegment const> l = std::dynamic_pointer_cast<LineSegment const>(c);

    if (l.get() == nullptr) {
        return false;
    } else {
        return (on_manifold(c->start()) && on_manifold(c->end()));
    }
}

double LineSegment::length() const { return hypot(end()->x() - start()->x(), end()->y() - start()->y()); };

void LineSegment::replace_verticies(std::vector<std::shared_ptr<Vertex const>> const &oldv, std::vector<std::shared_ptr<Vertex const>> const &newv) {
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