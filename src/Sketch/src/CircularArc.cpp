#include "Sketch.hpp"

void CircularArc::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    /*
        Calculate residual and jacobian from (1 - rad / radius())
        Normalize linearized equations by multiplying through by radius()
    */

    // #TODO: Clean this method
    double dx, dy, rad;

    // Jacobian, Residual
    dx = Start->x() - Center->x();
    dy = Start->y() - Center->y();
    rad = sqrt(dx * dx + dy * dy);

    r(EquationIndex) = radius() - rad;

    J(EquationIndex, Radius->get_index()) += rad / radius();

    J(EquationIndex, Start->X->get_index()) -= dx / rad;
    J(EquationIndex, Start->Y->get_index()) -= dy / rad;

    J(EquationIndex, Center->X->get_index()) += dx / rad;
    J(EquationIndex, Center->Y->get_index()) += dy / rad;

    dx = End->x() - Center->x();
    dy = End->y() - Center->y();
    rad = sqrt(dx * dx + dy * dy);

    r(EquationIndex + 1) = radius() - rad;

    J(EquationIndex + 1, Radius->get_index()) += rad / radius();

    J(EquationIndex + 1, End->X->get_index()) -= dx / rad;
    J(EquationIndex + 1, End->Y->get_index()) -= dy / rad;

    J(EquationIndex + 1, Center->X->get_index()) += dx / rad;
    J(EquationIndex + 1, Center->Y->get_index()) += dy / rad;
}

double CircularArc::s_to_a(double s) const {
    const double xc = Center->x();
    const double yc = Center->y();
    const double x0 = Start->x();
    const double y0 = Start->y();
    const double x1 = End->x();
    const double y1 = End->y();

    double a0 = atan2(y0 - yc, x0 - xc);
    double a1 = atan2(y1 - yc, x1 - xc);
    if (a1 < a0) {
        a1 += 2.0 * M_PI;
    }

    return (a0 * (1.0 - s) + a1 * s);
}

double CircularArc::a_to_s(double a) const {
    const double xc = Center->x();
    const double yc = Center->y();
    const double x0 = Start->x();
    const double y0 = Start->y();
    const double x1 = End->x();
    const double y1 = End->y();

    double a0 = atan2(y0 - yc, x0 - xc);
    double a1 = atan2(y1 - yc, x1 - xc);
    if (a1 < a0) {
        a1 += 2.0 * M_PI;
    }

    while (a < a0) {
        a += 2.0 * M_PI;
    }

    return (a - a0) / (a1 - a0);
}

double CircularArc::arc_angle() const {
    const double xc = Center->x();
    const double yc = Center->y();
    const double x0 = Start->x();
    const double y0 = Start->y();
    const double x1 = End->x();
    const double y1 = End->y();

    double a0 = atan2(y0 - yc, x0 - xc);
    double a1 = atan2(y1 - yc, x1 - xc);
    if (a1 < a0) {
        a1 += 2.0 * M_PI;
    }

    return (a1 - a0);
}

Vertex CircularArc::point(double s) const {
    double a = s_to_a(s);

    return Vertex{Center->x() + radius() * cos(a), Center->y() + radius() * sin(a)};
}

Vertex CircularArc::tangent(double s, bool orientation) const {
    double a = s_to_a(s);

    if (orientation) {
        return Vertex{-sin(a), cos(a)};
    } else {
        return Vertex{sin(a), -cos(a)};
    }
}

double CircularArc::area() const {
    double a = arc_angle();

    return 0.5 * (a - sin(a)) * radius() * radius();
}

double CircularArc::a(double s, bool orientation) const {
    double a = s_to_a(s);

    if (orientation) {
        if (a < M_PI_2) {
            a += M_PI_2;
        } else {
            a -= M_PI;
            a -= M_PI_2;
        }
    } else {
        if (a > -M_PI_2) {
            a -= M_PI_2;
        } else {
            a += M_PI;
            a += M_PI_2;
        }
    }

    return a;
}

double CircularArc::da(double s, bool orientation) const {
    // #TODO, Maybe should calculate this by taking the derivative of the angle in the gradient direction and projection, present implementation may cause errors in the future
    double sgn = SIGN(arc_angle());

    if (orientation) {
        return (sgn / radius());
    } else {
        return (-sgn / radius());
    }
}

double CircularArc::supremum() const {
    double sup = std::fmax(Start->hypot(), End->hypot());
    double s = a_to_s(Center->atan());

    if (s > 0 && s < 1) {
        sup = std::fmax(sup, point(s).hypot());
    }

    return std::move(sup);
}

bool CircularArc::on_manifold(const double x, const double y) const {
    double dx = Center->x() - x;
    double dy = Center->y() - y;
    double dr = sqrt(dx * dx + dy * dy);
    double tol = FLT_EPSILON * radius();

    if (abs(dr - radius()) < tol) {
        return true;
    } else {
        return false;
    }
}

bool CircularArc::is_identical(const Curve *c) const {
    const CircularArc *cc = dynamic_cast<const CircularArc *>(c);

    if (cc == nullptr) {
        return false;
    } else {
        return is_identical(
                cc->radius(),
                cc->center()->x(),
                cc->center()->y(),
                cc->start()->x(),
                cc->start()->y(),
                cc->end()->x(),
                cc->end()->y());
    }
}

bool CircularArc::is_identical(const Curve *c, const Vertex *origin, const double angle) const {
    const CircularArc *cc = dynamic_cast<const CircularArc *>(c);

    if (cc == nullptr) {
        return false;
    } else {
        double xc, yc, xs, ys, xe, ye;

        std::tie(xc, yc) = cc->center()->rotate(origin, angle);
        std::tie(xs, ys) = cc->start()->rotate(origin, angle);
        std::tie(xe, ye) = cc->end()->rotate(origin, angle);

        return is_identical(cc->radius(), xc, yc, xs, ys, xe, ye);
    }
}

bool CircularArc::is_identical(const double r, const double xc, const double yc, const double xs, const double ys,
                               const double xe, const double ye) const {
    double tol = FLT_EPSILON * radius();

    return abs(radius() - r) < tol
           && abs(center()->x() - xc) < tol
           && abs(center()->y() - yc) < tol
           && abs(start()->x() - xs) < tol
           && abs(start()->y() - ys) < tol
           && abs(end()->x() - xe) < tol
           && abs(end()->y() - ye) < tol;
}

bool CircularArc::is_coincident(const Curve *c) const {
    const CircularArc *cc = dynamic_cast<const CircularArc *>(c);

    if (cc == nullptr) {
        return false;
    } else {
        // #TODO: Extract input curve center, put the rest of the method in subroutine
        // #TODO: Then, can call subroutine for implementation of rotation version
        double xc = cc->Center->x();
        double yc = cc->Center->y();
        double rc = cc->radius();
        double tol = radius() * FLT_EPSILON;

        if (abs(Center->x() - xc) < tol && abs(Center->y() - yc) < tol && abs(radius() - rc) < tol) {
            return true;
        } else {
            return false;
        }
    }
}

bool CircularArc::on_segment(const double x, const double y) const {
    double dx = x - center()->x();
    double dy = y - center()->y();
    double dr = sqrt(dx * dx + dy * dy);
    double tol = radius() * FLT_EPSILON;

    if (abs(dr - radius()) < tol) {
        double da = atan2(dy, dx);

        double dx0 = start()->x() - center()->x();
        double dy0 = start()->y() - center()->y();
        double da0 = atan2(dy0, dx0);

        double dx1 = end()->x() - center()->x();
        double dy1 = end()->y() - center()->y();
        double da1 = atan2(dy1, dx1);

        if (da1 < da0) {
            da1 += 2.0 * M_PI;
        }

        if (da0 - da >= tol) {
            da += 2.0 * M_PI;
        }

        tol = FLT_EPSILON;
        if (da0 - da < tol && da - da1 < tol) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

void CircularArc::replace_verticies(std::vector<Vertex *> oldv, std::vector<Vertex *> newv) {
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

    i = std::find(oldv.begin(), oldv.end(), Center);
    if (i != oldv.end()) {
        size_t j = i - oldv.begin();
        Center = newv[j];
    }
}