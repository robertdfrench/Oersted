#include "Tangency.h"
#include "CircularArc.h"
#include "LineSegment.h"
#include "Vertex.h"

void Tangency::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
    /*
        Draw two vectors from the center of the circular arc to the end-points of the line
        segment. The line is tangent to the circle iff the cross product of the vectors is
        equal to the product of the length of the line and the radius of the arc, to within a
        sign. The cross product of the vectors is equal to twice the signed area of the
        triangle defined by the three points. The length of the line is equal to the base of
        the triangle. The radius of the circular is the height of the triangle.
    */

    double xc = Arc->center()->x();
    double yc = Arc->center()->y();
    double rc = Arc->radius();

    double x0 = Line->start()->x();
    double y0 = Line->start()->y();

    double x1 = Line->end()->x();
    double y1 = Line->end()->y();

    double dx = x1 - x0;
    double dy = y1 - y0;
    double dr = sqrt(dx * dx + dy * dy);

    double v0cx = x0 - xc;
    double v0cy = y0 - yc;

    double v1cx = x1 - xc;
    double v1cy = y1 - yc;

    double cross = (v0cx * v1cy - v0cy * v1cx) / dr;
    double f0, f1;

    dx /= dr;
    dy /= dr;
    dr *= SIGN(cross);

    r(EquationIndex) = abs(cross) - rc;

    J(EquationIndex, Arc->radius_index()) -= abs(cross) / rc;

    f0 = (v1cy + cross * dx) / dr;
    f1 = (v0cy + cross * dx) / dr;

    J(EquationIndex, Arc->center()->x_index()) -= f0 - f1;
    J(EquationIndex, Line->start()->x_index()) += f0;
    J(EquationIndex, Line->end()->x_index()) -= f1;

    f0 = (v1cx - cross * dy) / dr;
    f1 = (v0cx - cross * dy) / dr;

    J(EquationIndex, Arc->center()->y_index()) += f0 - f1;
    J(EquationIndex, Line->start()->y_index()) -= f0;
    J(EquationIndex, Line->end()->y_index()) += f1;
}