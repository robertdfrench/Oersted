#include "Distance.h"
#include "CircularArc.h"
#include "LineSegment.h"
#include "Vertex.h"

template<>
size_t Distance<CircularArc>::set_equation_index(size_t i) {
    EquationIndex = i;
    return 1;
};

template<>
void Distance<CircularArc>::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
    /*
        There are two different cases depending on whether or not the smaller circle is intended to be
        exterior or interior to the larger circle. At present, the intent is discriminated by the position
        of the smaller circle origin relative to the larger circle. The interior behavior is unstable in
        the sense that for large enough distances, the smaller circle will get pushed outside the larger
        circle and find a new equilibrium there.
    */

    double r0 = Element0->radius();
    double x0 = Element0->center()->x();
    double y0 = Element0->center()->y();

    double r1 = Element1->radius();
    double x1 = Element1->center()->x();
    double y1 = Element1->center()->y();

    double dx = x1 - x0;
    double dy = y1 - y0;
    double dr = sqrt(dx * dx + dy * dy);

    if (dr >= std::fmax(r0, r1)) {
        // Exterior
        r(EquationIndex) = dr - r0 - r1 - Dim;
        J(EquationIndex, Element0->radius_index()) -= (dr - Dim) / (r0 + r1);
        J(EquationIndex, Element1->radius_index()) -= (dr - Dim) / (r0 + r1);
    } else {
        // Interior
        if (r0 > r1) {
            r(EquationIndex) = Dim - (r0 - r1 - dr);
            J(EquationIndex, Element0->radius_index()) -= (Dim + dr) / (r0 - r1);
            J(EquationIndex, Element1->radius_index()) += (Dim + dr) / (r0 - r1);
        } else {
            r(EquationIndex) = Dim - (r1 - r0 - dr);
            J(EquationIndex, Element0->radius_index()) += (Dim + dr) / (r1 - r0);
            J(EquationIndex, Element1->radius_index()) -= (Dim + dr) / (r1 - r0);
        }
    }

    dx /= dr;
    J(EquationIndex, Element0->center()->x_index()) -= dx;
    J(EquationIndex, Element1->center()->x_index()) += dx;

    dy /= dr;
    J(EquationIndex, Element0->center()->y_index()) -= dy;
    J(EquationIndex, Element1->center()->y_index()) += dy;
}

template
class Distance<CircularArc>;

template<>
size_t Distance<LineSegment>::set_equation_index(size_t i) {
    EquationIndex = i;
    return 2;
};

template<>
void Distance<LineSegment>::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
    /*
        Pick one line, draw vectors from the center of the line to each of the remaining line's end points.
        The area of the triangle formed by the vectors and line should equal half of the length of the line times
        the distance beween the lines (one half base times height). The signed area of the triangle is equal to
        the cross product of the two vectors.
    */
    double x00 = Element0->start()->x();
    double y00 = Element0->start()->y();
    double x01 = Element0->end()->x();
    double y01 = Element0->end()->y();
    double x0m = 0.5 * (x00 + x01);
    double y0m = 0.5 * (y00 + y01);

    double v0x = x01 - x00;
    double v0y = y01 - y00;
    double d0 = sqrt(v0x * v0x + v0y * v0y);
    v0x /= d0;
    v0y /= d0;

    double x10 = Element1->start()->x();
    double y10 = Element1->start()->y();
    double x11 = Element1->end()->x();
    double y11 = Element1->end()->y();
    double x1m = 0.5 * (x10 + x11);
    double y1m = 0.5 * (y10 + y11);

    double v1x = x11 - x10;
    double v1y = y11 - y10;
    double d1 = sqrt(v1x * v1x + v1y * v1y);
    v1x /= d1;
    v1y /= d1;

    double v0mx, v0my, v1mx, v1my, cross, dot, d01, scale, f;

    // Element0 center to Element 1 end points
    v0mx = x10 - x0m;
    v0my = y10 - y0m;
    v1mx = x11 - x0m;
    v1my = y11 - y0m;

    d01 = d0 + d1;
    cross = (v0mx * v1my - v1mx * v0my) / d01;
    d01 *= SIGN(cross);

    r(EquationIndex) = abs(cross) - Dim / 2.0;

    J(EquationIndex, Element0->start()->x_index()) -= (v1y * d1 - cross * v0x) / d01;
    J(EquationIndex, Element0->end()->x_index()) -= (v1y * d1 + cross * v0x) / d01;

    J(EquationIndex, Element1->start()->x_index()) += (v1my + cross * v1x) / d01;
    J(EquationIndex, Element1->end()->x_index()) -= (v0my + cross * v1x) / d01;

    J(EquationIndex, Element0->start()->y_index()) += (v1x * d1 + cross * v0y) / d01;
    J(EquationIndex, Element0->end()->y_index()) += (v1x * d1 - cross * v0y) / d01;

    J(EquationIndex, Element1->start()->y_index()) -= (v1mx - cross * v1y) / d01;
    J(EquationIndex, Element1->end()->y_index()) += (v0mx - cross * v1y) / d01;

    // Element1 center to Element 0 end points
    v0mx = x00 - x1m;
    v0my = y00 - y1m;
    v1mx = x01 - x1m;
    v1my = y01 - y1m;

    d01 = d0 + d1;
    cross = (v0mx * v1my - v1mx * v0my) / d01;
    d01 *= SIGN(cross);

    r(EquationIndex) += abs(cross) - Dim / 2.0;

    J(EquationIndex, Element0->start()->x_index()) += (v1my + cross * v0x) / d01;
    J(EquationIndex, Element0->end()->x_index()) -= (v0my + cross * v0x) / d01;

    J(EquationIndex, Element1->start()->x_index()) -= (v0y * d0 - cross * v1x) / d01;
    J(EquationIndex, Element1->end()->x_index()) -= (v0y * d0 + cross * v1x) / d01;

    J(EquationIndex, Element0->start()->y_index()) -= (v1mx - cross * v0y) / d01;
    J(EquationIndex, Element0->end()->y_index()) += (v0mx - cross * v0y) / d01;

    J(EquationIndex, Element1->start()->y_index()) += (v0x * d0 + cross * v1y) / d01;
    J(EquationIndex, Element1->end()->y_index()) += (v0x * d0 - cross * v1y) / d01;

    // Parallel Constraint
    cross = (v0x * v1y - v0y * v1x);
    dot = (v0x * v1x + v0y * v1y);
    scale = std::fmax(d0, d1);

    if (abs(cross) < abs(dot)) {
        // Use cross product equation
        r(EquationIndex + 1) = scale * cross;

        d0 /= scale;
        d1 /= scale;

        f = (v1y - cross * v0x) / d0;
        J(EquationIndex + 1, Element0->start()->x_index()) -= f;
        J(EquationIndex + 1, Element0->end()->x_index()) += f;

        f = (v1x + cross * v0y) / d0;
        J(EquationIndex + 1, Element0->start()->y_index()) += f;
        J(EquationIndex + 1, Element0->end()->y_index()) -= f;

        f = (v0y + cross * v1x) / d1;
        J(EquationIndex + 1, Element1->start()->x_index()) += f;
        J(EquationIndex + 1, Element1->end()->x_index()) -= f;

        f = (v0x - cross * v1y) / d1;
        J(EquationIndex + 1, Element1->start()->y_index()) -= f;
        J(EquationIndex + 1, Element1->end()->y_index()) += f;
    } else {
        // Use dot product equation
        r(EquationIndex + 1) = scale * (abs(dot) - 1.0);

        d0 /= (scale * SIGN(dot));
        d1 /= (scale * SIGN(dot));

        f = (v1x - dot * v0x) / d0;
        J(EquationIndex + 1, Element0->start()->x_index()) -= f;
        J(EquationIndex + 1, Element0->end()->x_index()) += f;

        f = (v1y - dot * v0y) / d0;
        J(EquationIndex + 1, Element0->start()->y_index()) -= f;
        J(EquationIndex + 1, Element0->end()->y_index()) += f;

        f = (v0x - dot * v1x) / d1;
        J(EquationIndex + 1, Element1->start()->x_index()) -= f;
        J(EquationIndex + 1, Element1->end()->x_index()) += f;

        f = (v0y - dot * v1y) / d1;
        J(EquationIndex + 1, Element1->start()->y_index()) -= f;
        J(EquationIndex + 1, Element1->end()->y_index()) += f;
    }
}

template
class Distance<LineSegment>;

template<>
size_t Distance<Vertex>::set_equation_index(size_t i) {
    EquationIndex = i;
    return 1;
}

template<>
void Distance<Vertex>::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
    double x0 = Element0->x();
    double y0 = Element0->y();
    double x1 = Element1->x();
    double y1 = Element1->y();

    double dx = x1 - x0;
    double dy = y1 - y0;
    double dr = sqrt(dx * dx + dy * dy);
    dx = dx / dr;
    dy = dy / dr;

    r(EquationIndex) = dr - Dim;

    J(EquationIndex, Element0->x_index()) -= dx;
    J(EquationIndex, Element0->y_index()) -= dy;
    J(EquationIndex, Element1->x_index()) += dx;
    J(EquationIndex, Element1->y_index()) += dy;
}

template
class Distance<Vertex>;