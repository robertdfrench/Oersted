#include "Sketch.hpp"

// #TODO: Split into multiple cpp files
// Angle
void Angle::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    /*
        The dot product of two unit vectors equals the cosine of the angle between them.
        The cross product of two unit vectros equals the sine of the angle between them.
        To solve the constraint, use the smaller of the dot produce and cross product equations.
        Do this because the derivative of the dot (resp. cross) product with respect to the
        parameters is zero when the cosine (resp. sine) is equal to one.
    */

    const double x00 = Line0->Start->x();
    const double x01 = Line0->End->x();
    const double x10 = Line1->Start->x();
    const double x11 = Line1->End->x();

    const double y00 = Line0->Start->y();
    const double y01 = Line0->End->y();
    const double y10 = Line1->Start->y();
    const double y11 = Line1->End->y();

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
        J(EquationIndex, Line0->Start->X->get_index()) -= f;
        J(EquationIndex, Line0->End->X->get_index()) += f;

        f = scale * (vx0 - dot * vx1) / d1;
        J(EquationIndex, Line1->Start->X->get_index()) -= f;
        J(EquationIndex, Line1->End->X->get_index()) += f;

        f = scale * (vy1 - dot * vy0) / d0;
        J(EquationIndex, Line0->Start->Y->get_index()) -= f;
        J(EquationIndex, Line0->End->Y->get_index()) += f;

        f = scale * (vy0 - dot * vy1) / d1;
        J(EquationIndex, Line1->Start->Y->get_index()) -= f;
        J(EquationIndex, Line1->End->Y->get_index()) += f;
    } else {
        // Solve (cross product) = sin(Dim)
        r(EquationIndex) = scale * (cross - sin(M_PI * Dim / 180.0));

        f = scale * (vy1 - cross * vx0) / d0;
        J(EquationIndex, Line0->Start->X->get_index()) -= f;
        J(EquationIndex, Line0->End->X->get_index()) += f;

        f = scale * (vy0 + cross * vx1) / d1;
        J(EquationIndex, Line1->Start->X->get_index()) += f;
        J(EquationIndex, Line1->End->X->get_index()) -= f;

        f = scale * (vx1 + cross * vy0) / d0;
        J(EquationIndex, Line0->Start->Y->get_index()) += f;
        J(EquationIndex, Line0->End->Y->get_index()) -= f;

        f = scale * (vx0 - cross * vy1) / d1;
        J(EquationIndex, Line1->Start->Y->get_index()) -= f;
        J(EquationIndex, Line1->End->Y->get_index()) += f;
    }
}

// Coincident
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

// Distance
template<>
size_t Distance<CircularArc>::set_equation_index(size_t i) {
    EquationIndex = i;
    return 1;
};

template<>
void Distance<CircularArc>::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    /*
        There are two different cases depending on whether or not the smaller circle is intended to be
        exterior or interior to the larger circle. At present, the intent is discriminated by the position
        of the smaller circle origin relative to the larger circle. The interior behavior is unstable in
        the sense that for large enough distances, the smaller circle will get pushed outside the larger
        circle and find a new equilibrium there.
    */

    const double r0 = Element0->radius();
    const double x0 = Element0->Center->x();
    const double y0 = Element0->Center->y();

    const double r1 = Element1->radius();
    const double x1 = Element1->Center->x();
    const double y1 = Element1->Center->y();

    double dx = x1 - x0;
    double dy = y1 - y0;
    double dr = sqrt(dx * dx + dy * dy);

    if (dr >= std::fmax(r0, r1)) {
        // Exterior
        r(EquationIndex) = dr - r0 - r1 - Dim;
        J(EquationIndex, Element0->Radius->get_index()) -= (dr - Dim) / (r0 + r1);
        J(EquationIndex, Element1->Radius->get_index()) -= (dr - Dim) / (r0 + r1);
    } else {
        // Interior
        if (r0 > r1) {
            r(EquationIndex) = Dim - (r0 - r1 - dr);
            J(EquationIndex, Element0->Radius->get_index()) -= (Dim + dr) / (r0 - r1);
            J(EquationIndex, Element1->Radius->get_index()) += (Dim + dr) / (r0 - r1);
        } else {
            r(EquationIndex) = Dim - (r1 - r0 - dr);
            J(EquationIndex, Element0->Radius->get_index()) += (Dim + dr) / (r1 - r0);
            J(EquationIndex, Element1->Radius->get_index()) -= (Dim + dr) / (r1 - r0);
        }
    }

    dx /= dr;
    J(EquationIndex, Element0->Center->X->get_index()) -= dx;
    J(EquationIndex, Element1->Center->X->get_index()) += dx;

    dy /= dr;
    J(EquationIndex, Element0->Center->Y->get_index()) -= dy;
    J(EquationIndex, Element1->Center->Y->get_index()) += dy;
}

template
class Distance<CircularArc>;

template<>
size_t Distance<LineSegment>::set_equation_index(size_t i) {
    EquationIndex = i;
    return 2;
};

template<>
void Distance<LineSegment>::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    /*
        Pick one line, draw vectors from the center of the line to each of the remaining line's end points.
        The area of the triangle formed by the vectors and line should equal half of the length of the line times
        the distance beween the lines (one half base times height). The signed area of the triangle is equal to
        the cross product of the two vectors.
    */
    const double x00 = Element0->Start->x();
    const double y00 = Element0->Start->y();
    const double x01 = Element0->End->x();
    const double y01 = Element0->End->y();
    const double x0m = 0.5 * (x00 + x01);
    const double y0m = 0.5 * (y00 + y01);

    double v0x = x01 - x00;
    double v0y = y01 - y00;
    double d0 = sqrt(v0x * v0x + v0y * v0y);
    v0x /= d0;
    v0y /= d0;

    const double x10 = Element1->Start->x();
    const double y10 = Element1->Start->y();
    const double x11 = Element1->End->x();
    const double y11 = Element1->End->y();
    const double x1m = 0.5 * (x10 + x11);
    const double y1m = 0.5 * (y10 + y11);

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

    J(EquationIndex, Element0->Start->X->get_index()) -= (v1y * d1 - cross * v0x) / d01;
    J(EquationIndex, Element0->End->X->get_index()) -= (v1y * d1 + cross * v0x) / d01;

    J(EquationIndex, Element1->Start->X->get_index()) += (v1my + cross * v1x) / d01;
    J(EquationIndex, Element1->End->X->get_index()) -= (v0my + cross * v1x) / d01;

    J(EquationIndex, Element0->Start->Y->get_index()) += (v1x * d1 + cross * v0y) / d01;
    J(EquationIndex, Element0->End->Y->get_index()) += (v1x * d1 - cross * v0y) / d01;

    J(EquationIndex, Element1->Start->Y->get_index()) -= (v1mx - cross * v1y) / d01;
    J(EquationIndex, Element1->End->Y->get_index()) += (v0mx - cross * v1y) / d01;

    // Element1 center to Element 0 end points
    v0mx = x00 - x1m;
    v0my = y00 - y1m;
    v1mx = x01 - x1m;
    v1my = y01 - y1m;

    d01 = d0 + d1;
    cross = (v0mx * v1my - v1mx * v0my) / d01;
    d01 *= SIGN(cross);

    r(EquationIndex) += abs(cross) - Dim / 2.0;

    J(EquationIndex, Element0->Start->X->get_index()) += (v1my + cross * v0x) / d01;
    J(EquationIndex, Element0->End->X->get_index()) -= (v0my + cross * v0x) / d01;

    J(EquationIndex, Element1->Start->X->get_index()) -= (v0y * d0 - cross * v1x) / d01;
    J(EquationIndex, Element1->End->X->get_index()) -= (v0y * d0 + cross * v1x) / d01;

    J(EquationIndex, Element0->Start->Y->get_index()) -= (v1mx - cross * v0y) / d01;
    J(EquationIndex, Element0->End->Y->get_index()) += (v0mx - cross * v0y) / d01;

    J(EquationIndex, Element1->Start->Y->get_index()) += (v0x * d0 + cross * v1y) / d01;
    J(EquationIndex, Element1->End->Y->get_index()) += (v0x * d0 - cross * v1y) / d01;

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
        J(EquationIndex + 1, Element0->Start->X->get_index()) -= f;
        J(EquationIndex + 1, Element0->End->X->get_index()) += f;

        f = (v1x + cross * v0y) / d0;
        J(EquationIndex + 1, Element0->Start->Y->get_index()) += f;
        J(EquationIndex + 1, Element0->End->Y->get_index()) -= f;

        f = (v0y + cross * v1x) / d1;
        J(EquationIndex + 1, Element1->Start->X->get_index()) += f;
        J(EquationIndex + 1, Element1->End->X->get_index()) -= f;

        f = (v0x - cross * v1y) / d1;
        J(EquationIndex + 1, Element1->Start->Y->get_index()) -= f;
        J(EquationIndex + 1, Element1->End->Y->get_index()) += f;
    } else {
        // Use dot product equation
        r(EquationIndex + 1) = scale * (abs(dot) - 1.0);

        d0 /= (scale * SIGN(dot));
        d1 /= (scale * SIGN(dot));

        f = (v1x - dot * v0x) / d0;
        J(EquationIndex + 1, Element0->Start->X->get_index()) -= f;
        J(EquationIndex + 1, Element0->End->X->get_index()) += f;

        f = (v1y - dot * v0y) / d0;
        J(EquationIndex + 1, Element0->Start->Y->get_index()) -= f;
        J(EquationIndex + 1, Element0->End->Y->get_index()) += f;

        f = (v0x - dot * v1x) / d1;
        J(EquationIndex + 1, Element1->Start->X->get_index()) -= f;
        J(EquationIndex + 1, Element1->End->X->get_index()) += f;

        f = (v0y - dot * v1y) / d1;
        J(EquationIndex + 1, Element1->Start->Y->get_index()) -= f;
        J(EquationIndex + 1, Element1->End->Y->get_index()) += f;
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
void Distance<Vertex>::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    const double x0 = Element0->x();
    const double y0 = Element0->y();
    const double x1 = Element1->x();
    const double y1 = Element1->y();

    double dx = x1 - x0;
    double dy = y1 - y0;
    double dr = sqrt(dx * dx + dy * dy);
    dx = dx / dr;
    dy = dy / dr;

    r(EquationIndex) = dr - Dim;

    J(EquationIndex, Element0->X->get_index()) -= dx;
    J(EquationIndex, Element0->Y->get_index()) -= dy;
    J(EquationIndex, Element1->X->get_index()) += dx;
    J(EquationIndex, Element1->Y->get_index()) += dy;
}

template
class Distance<Vertex>;

// Fixation
Fixation::Fixation(Vertex &v) {
    Point = &v;
    Dim = new Vertex(v);
}

void Fixation::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    r(EquationIndex) = Point->x() - Dim->x();
    J(EquationIndex, Point->X->get_index()) += 1.0;

    r(EquationIndex + 1) = Point->y() - Dim->y();
    J(EquationIndex + 1, Point->Y->get_index()) += 1.0;
}

// Horizontal
void Horizontal::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    r(EquationIndex) = Line->Start->y() - Line->End->y();

    J(EquationIndex, Line->Start->Y->get_index()) += 1.0;
    J(EquationIndex, Line->End->Y->get_index()) -= 1.0;
}

// Length
void Length::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    double dx = Line->End->x() - Line->Start->x();
    double dy = Line->End->y() - Line->Start->y();
    const double d = sqrt(dx * dx + dy * dy);

    r(EquationIndex) = d - Dim;

    dx = dx / d;
    dy = dy / d;

    J(EquationIndex, Line->Start->X->get_index()) -= dx;
    J(EquationIndex, Line->Start->Y->get_index()) -= dy;

    J(EquationIndex, Line->End->X->get_index()) += dx;
    J(EquationIndex, Line->End->Y->get_index()) += dy;
}

// Radius
void Radius::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    r(EquationIndex) = Arc->radius() - Dim;

    J(EquationIndex, Arc->Radius->get_index()) += 1.0;
}

// Rotation
void Rotation::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    double xo = Origin->x();
    double yo = Origin->y();

    double x0 = V0->x();
    double y0 = V0->y();

    double x1 = V1->x();
    double y1 = V1->y();

    double a = Angle * M_PI / 180.0;

    size_t ei = EquationIndex;
    r(ei) = (x0 - xo) * cos(a) - (y0 - yo) * sin(a) - (x1 - xo);

    J(ei, Origin->X->get_index()) = 1.0 - cos(a);
    J(ei, Origin->Y->get_index()) = sin(a);
    J(ei, V0->X->get_index()) = cos(a);
    J(ei, V0->Y->get_index()) = -sin(a);
    J(ei, V1->X->get_index()) = -1.0;

    ei++;
    r(ei) = (x0 - xo) * sin(a) + (y0 - yo) * cos(a) - (y1 - yo);

    J(ei, Origin->X->get_index()) = -sin(a);
    J(ei, Origin->Y->get_index()) = 1 - cos(a);
    J(ei, V0->X->get_index()) = sin(a);
    J(ei, V0->Y->get_index()) = cos(a);
    J(ei, V1->Y->get_index()) = -1.0;
}

// Symmetry
void Symmetry::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    double x0 = V0->x();
    double y0 = V0->y();

    double x1 = V1->x();
    double y1 = V1->y();

    double dx01 = x0 - x1;
    double dy01 = y0 - y1;

    double d01 = sqrt(dx01 * dx01 + dy01 * dy01);

    double xs = SymmetryLine->start()->x();
    double ys = SymmetryLine->start()->y();

    double xe = SymmetryLine->end()->x();
    double ye = SymmetryLine->end()->y();

    double dxes = xe - xs;
    double dyes = ye - ys;

    double dse = sqrt(dxes * dxes + dyes * dyes);

    double scale = 1.0 / sqrt(d01 * dse);

    r(EquationIndex) = ((x0 - x1) * (xe - xs) + (y0 - y1) * (ye - ys)) * scale;

    J(EquationIndex, SymmetryLine->start()->X->get_index()) -= dx01 * scale;
    J(EquationIndex, SymmetryLine->start()->Y->get_index()) -= dy01 * scale;

    J(EquationIndex, SymmetryLine->end()->X->get_index()) += dx01 * scale;
    J(EquationIndex, SymmetryLine->end()->Y->get_index()) += dy01 * scale;

    J(EquationIndex, V0->X->get_index()) += dxes * scale;
    J(EquationIndex, V0->Y->get_index()) += dyes * scale;

    J(EquationIndex, V1->X->get_index()) -= dxes * scale;
    J(EquationIndex, V1->Y->get_index()) -= dyes * scale;

    r(EquationIndex + 1) =
            (((x0 - xs) * (y0 - ye) - (x0 - xe) * (y0 - ys)) + ((x1 - xs) * (y1 - ye) - (x1 - xe) * (y1 - ys))) * scale;

    J(EquationIndex + 1, SymmetryLine->start()->X->get_index()) -= (y1 + y0 - 2.0 * ye) * scale;
    J(EquationIndex + 1, SymmetryLine->start()->Y->get_index()) += (x0 + x1 - 2.0 * xe) * scale;

    J(EquationIndex + 1, SymmetryLine->end()->X->get_index()) += (y1 + y0 - 2.0 * ye) * scale;
    J(EquationIndex + 1, SymmetryLine->end()->Y->get_index()) -= (x0 + x1 - 2.0 * xe) * scale;

    J(EquationIndex + 1, V0->X->get_index()) -= dyes * scale;
    J(EquationIndex + 1, V0->Y->get_index()) += dxes * scale;

    J(EquationIndex + 1, V1->X->get_index()) -= dyes * scale;
    J(EquationIndex + 1, V1->Y->get_index()) += dxes * scale;
}

// Tagency
void Tangency::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    /*
        Draw two vectors from the center of the circlular arc to the end-points of the line
        segment. The line is tagent to the circle iff the cross product of the vectors is
        equal to the product of the length of the line and the radius of the arc, to within a
        sign. The cross product of the vectors is equal to twice the signed area of the
        triangle defined by the three points. The length of the line is equal to the base of
        the triangle. The radius of the circular is the height of the triangle.
    */

    const double xc = Arc->Center->x();
    const double yc = Arc->Center->y();
    const double rc = Arc->radius();

    const double x0 = Line->Start->x();
    const double y0 = Line->Start->y();

    const double x1 = Line->End->x();
    const double y1 = Line->End->y();

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

    J(EquationIndex, Arc->Radius->get_index()) -= abs(cross) / rc;

    f0 = (v1cy + cross * dx) / dr;
    f1 = (v0cy + cross * dx) / dr;

    J(EquationIndex, Arc->Center->X->get_index()) -= f0 - f1;
    J(EquationIndex, Line->Start->X->get_index()) += f0;
    J(EquationIndex, Line->End->X->get_index()) -= f1;

    f0 = (v1cx - cross * dy) / dr;
    f1 = (v0cx - cross * dy) / dr;

    J(EquationIndex, Arc->Center->Y->get_index()) += f0 - f1;
    J(EquationIndex, Line->Start->Y->get_index()) -= f0;
    J(EquationIndex, Line->End->Y->get_index()) += f1;
}

// Vertical
void Vertical::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    r(EquationIndex) = Line->Start->x() - Line->End->x();

    J(EquationIndex, Line->Start->X->get_index()) += 1.0;
    J(EquationIndex, Line->End->X->get_index()) -= 1.0;
}