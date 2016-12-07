#include "Angle.h"
#include "LineSegment.h"

void Angle::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
    /*
        The dot product of two unit vectors equals the cosine of the angle between them.
        The cross product of two unit vectros equals the sine of the angle between them.
        To solve the constraint, use the smaller of the dot produce and cross product equations.
        Do this because the derivative of the dot (resp. cross) product with respect to the
        parameters is zero when the cosine (resp. sine) is equal to one.
    */

    double x00 = Line0->start()->x();
    double x01 = Line0->end()->x();
    double x10 = Line1->start()->x();
    double x11 = Line1->end()->x();

    double y00 = Line0->start()->y();
    double y01 = Line0->end()->y();
    double y10 = Line1->start()->y();
    double y11 = Line1->end()->y();

    double vx0 = x01 - x00;
    double vy0 = y01 - y00;
    double vx1 = x11 - x10;
    double vy1 = y11 - y10;

    double d0 = sqrt(vx0 * vx0 + vy0 * vy0);
    double d1 = sqrt(vx1 * vx1 + vy1 * vy1);

    vx0 /= d0;
    vy0 /= d0;
    vx1 /= d1;
    vy1 /= d1;

    double dot = (vx0 * vx1 + vy0 * vy1);
    double cross = (vx0 * vy1 - vy0 * vx1);
    double scale = std::fmax(d0, d1);
    double f;

    if (abs(dot) < abs(cross)) {
        // Solve (dot product) = cos(Dim)
        r(EquationIndex) = scale * (dot - cos(M_PI * Dim / 180.0));

        f = scale * (vx1 - dot * vx0) / d0;
        J(EquationIndex, Line0->start()->x_index()) -= f;
        J(EquationIndex, Line0->end()->x_index()) += f;

        f = scale * (vx0 - dot * vx1) / d1;
        J(EquationIndex, Line1->start()->x_index()) -= f;
        J(EquationIndex, Line1->end()->x_index()) += f;

        f = scale * (vy1 - dot * vy0) / d0;
        J(EquationIndex, Line0->start()->y_index()) -= f;
        J(EquationIndex, Line0->end()->y_index()) += f;

        f = scale * (vy0 - dot * vy1) / d1;
        J(EquationIndex, Line1->start()->y_index()) -= f;
        J(EquationIndex, Line1->end()->y_index()) += f;
    } else {
        // Solve (cross product) = sin(Dim)
        r(EquationIndex) = scale * (cross - sin(M_PI * Dim / 180.0));

        f = scale * (vy1 - cross * vx0) / d0;
        J(EquationIndex, Line0->start()->x_index()) -= f;
        J(EquationIndex, Line0->end()->x_index()) += f;

        f = scale * (vy0 + cross * vx1) / d1;
        J(EquationIndex, Line1->start()->x_index()) += f;
        J(EquationIndex, Line1->end()->x_index()) -= f;

        f = scale * (vx1 + cross * vy0) / d0;
        J(EquationIndex, Line0->start()->y_index()) += f;
        J(EquationIndex, Line0->end()->y_index()) -= f;

        f = scale * (vx0 - cross * vy1) / d1;
        J(EquationIndex, Line1->start()->y_index()) -= f;
        J(EquationIndex, Line1->end()->y_index()) += f;
    }
}