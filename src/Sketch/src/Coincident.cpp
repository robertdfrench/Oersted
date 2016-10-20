#include "Sketch.hpp"

template<>
void Coincident<CircularArc>::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    const double rc = Element->radius();
    const double xc = Element->Center->x();
    const double yc = Element->Center->y();
    const double xp = Point->x();
    const double yp = Point->y();

    double dx = xp - xc;
    double dy = yp - yc;
    double dr = sqrt(dx * dx + dy * dy);

    dx /= dr;
    dy /= dr;

    r(EquationIndex) = dr - rc;

    J(EquationIndex, Element->Radius->get_index()) = -dr / rc;

    J(EquationIndex, Element->Center->X->get_index()) -= dx;
    J(EquationIndex, Element->Center->Y->get_index()) -= dy;

    J(EquationIndex, Point->X->get_index()) += dx;
    J(EquationIndex, Point->Y->get_index()) += dy;
}

template
class Coincident<CircularArc>;

template<>
void Coincident<LineSegment>::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    const double xp = Point->x();
    const double yp = Point->y();

    const double x0 = Element->Start->x();
    const double y0 = Element->Start->y();
    const double x1 = Element->End->x();
    const double y1 = Element->End->y();

    double dx0 = x0 - xp;
    double dy0 = y0 - yp;
    double dr0 = sqrt(dx0 * dx0 + dy0 * dy0);

    double dx1 = x1 - xp;
    double dy1 = y1 - yp;
    double dr1 = sqrt(dx1 * dx1 + dy1 * dy1);

    dx0 /= dr0;
    dy0 /= dr0;
    dx1 /= dr1;
    dy1 /= dr1;

    double cross = (dx0 * dy1 - dy0 * dx1);
    double scale = std::fmax(dr0, dr1);

    r(EquationIndex) = scale * cross;

    double f0, f1;

    f0 = scale * (dy1 - cross * dx0) / dr0;
    f1 = scale * (dy0 + cross * dx1) / dr1;

    J(EquationIndex, Point->X->get_index()) -= f0 - f1;
    J(EquationIndex, Element->Start->X->get_index()) += f0;
    J(EquationIndex, Element->End->X->get_index()) -= f1;

    f0 = scale * (dx1 + cross * dy0) / dr0;
    f1 = scale * (dx0 - cross * dy1) / dr1;

    J(EquationIndex, Point->Y->get_index()) += f0 - f1;
    J(EquationIndex, Element->Start->Y->get_index()) -= f0;
    J(EquationIndex, Element->End->Y->get_index()) += f1;
}

template
class Coincident<LineSegment>;