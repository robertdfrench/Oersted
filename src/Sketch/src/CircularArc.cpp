#include "CircularArc.h"
#include "doublen.h"

void CircularArc::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
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

    J(EquationIndex, radius_index()) += rad / radius();

    J(EquationIndex, Start->x_index()) -= dx / rad;
    J(EquationIndex, Start->y_index()) -= dy / rad;

    J(EquationIndex, Center->x_index()) += dx / rad;
    J(EquationIndex, Center->y_index()) += dy / rad;

    dx = End->x() - Center->x();
    dy = End->y() - Center->y();
    rad = sqrt(dx * dx + dy * dy);

    r(EquationIndex + 1) = radius() - rad;

    J(EquationIndex + 1, radius_index()) += rad / radius();

    J(EquationIndex + 1, End->x_index()) -= dx / rad;
    J(EquationIndex + 1, End->y_index()) -= dy / rad;

    J(EquationIndex + 1, Center->x_index()) += dx / rad;
    J(EquationIndex + 1, Center->y_index()) += dy / rad;
}

double CircularArc::s_to_a(double s) const {
    double xc = Center->x();
    double yc = Center->y();
    double x0 = Start->x();
    double y0 = Start->y();
    double x1 = End->x();
    double y1 = End->y();

    double a0 = atan2(y0 - yc, x0 - xc);
    double a1 = atan2(y1 - yc, x1 - xc);
    if (a1 < a0) {
        a1 += 2.0 * M_PI;
    }

    return (a0 * (1.0 - s) + a1 * s);
}

double CircularArc::a_to_s(double a) const {
    double xc = Center->x();
    double yc = Center->y();
    double x0 = Start->x();
    double y0 = Start->y();
    double x1 = End->x();
    double y1 = End->y();

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
    double xc = Center->x();
    double yc = Center->y();
    double x0 = Start->x();
    double y0 = Start->y();
    double x1 = End->x();
    double y1 = End->y();

    double a0 = atan2(y0 - yc, x0 - xc);
    double a1 = atan2(y1 - yc, x1 - xc);
    if (a1 < a0) {
        a1 += 2.0 * M_PI;
    }

    return (a1 - a0);
}

double2 CircularArc::point(double s) const {
    double a = s_to_a(s);

    return double2{center()->x() + radius() * cos(a), center()->y() + radius() * sin(a)};
}

double2 CircularArc::tangent(double s, bool orientation) const {
    double a = s_to_a(s);

    if (orientation) {
        return double2{-sin(a), cos(a)};
    } else {
        return double2{sin(a), -cos(a)};
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
    // TODO: Maybe should calculate this by taking the derivative of the angle in the gradient direction and projection, present implementation may cause errors in the future
    double sgn = SIGN(arc_angle());

    if (orientation) {
        return (sgn / radius());
    } else {
        return (-sgn / radius());
    }
}

double2 CircularArc::supremum() const {
    double x = start()->x();
    double y = start()->y();
    double sup = sqrt(x * x + y * y);
    double par = 0.0;

    double xx = end()->x();
    double yy = end()->y();
    double val = sqrt(xx * xx + yy * yy);
    if (val > sup) {
        x = xx;
        y = yy;
        sup = val;
        par = 1.0;
    }

    double s = a_to_s(center()->atan());
    if (s > 0.0 && s < 1.0) {
        double2 v = point(s);
        xx = v.X;
        yy = v.Y;
        val = sqrt(xx * xx + yy * yy);

        if (val > sup) {
            x = xx;
            y = yy;
            sup = val;
            par = s;
        }
    }

    double ang = s_to_a(par);
    double cross = abs(x * cos(ang) + y * sin(ang)) / sup; // cross product of vector from origin to point and tangent vector

    return double2{sup, cross};
}

bool CircularArc::on_manifold(double x, double y) const {
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

MatchOrientation CircularArc::is_identical(std::shared_ptr<Curve const> const &c) const {
    std::shared_ptr<CircularArc const> cc = std::dynamic_pointer_cast<CircularArc const>(c);

    if (cc.get() == nullptr) {
        return MatchOrientation::None;
    } else {
        return is_identical(cc->radius(), cc->center()->x(), cc->center()->y(), cc->start()->x(), cc->start()->y(), cc->end()->x(), cc->end()->y());
    }
}

MatchOrientation CircularArc::is_identical(std::shared_ptr<Curve const> const &c, std::shared_ptr<Vertex const> const &origin, double angle) const {
    std::shared_ptr<CircularArc const> cc = std::dynamic_pointer_cast<CircularArc const>(c);

    if (cc.get() == nullptr) {
        return MatchOrientation::None;
    } else {
        double2 pc = cc->center()->rotate(origin, angle);
        double2 ps = cc->start()->rotate(origin, angle);
        double2 pe = cc->end()->rotate(origin, angle);

        return is_identical(cc->radius(), pc.X, pc.Y, ps.X, ps.Y, pe.X, pe.Y);
    }
}

MatchOrientation CircularArc::is_identical(double r, double xc, double yc, double xs, double ys, double xe, double ye) const {
    double tol = FLT_EPSILON * radius();

    if (abs(radius() - r) < tol && abs(center()->x() - xc) < tol && abs(center()->y() - yc) < tol &&
        abs(start()->x() - xs) < tol && abs(start()->y() - ys) < tol && abs(end()->x() - xe) < tol &&
        abs(end()->y() - ye) < tol) {
        return MatchOrientation::Forward;
    } else {
        return MatchOrientation::None;
    }
}

bool CircularArc::is_coincident(std::shared_ptr<Curve const> const &c) const {
    std::shared_ptr<CircularArc const> cc = std::dynamic_pointer_cast<CircularArc const>(c);

    if (cc.get() == nullptr) {
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

bool CircularArc::on_segment(double x, double y) const {
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

void CircularArc::replace_verticies(std::vector<std::shared_ptr<Vertex const>> const &oldv, std::vector<std::shared_ptr<Vertex const>> const &newv) {
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