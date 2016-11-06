#include "Mesh.hpp"

// TODO: Make these mesh routines
bool are_intersecting(Edge const *e0, Edge const *e1, Mesh const &m) {
    // TODO, Make more detailed return type enumeration
    if (e0->ConstraintCurve != nullptr && e0->ConstraintCurve == e1->ConstraintCurve) {
        return false;
    }

    Point const v00 = m.point(e0->base());
    Point const v01 = m.point(e0->tip());
    Point const v10 = m.point(e1->base());
    Point const v11 = m.point(e1->tip());

    double xs0 = (v00.X + v01.X) / 2.0;
    double ys0 = (v00.Y + v01.Y) / 2.0;
    double xs1 = (v10.X + v11.X) / 2.0;
    double ys1 = (v10.Y + v11.Y) / 2.0;

    double xd0 = (v00.X - v01.X) / 2.0;
    double yd0 = (v00.Y - v01.Y) / 2.0;
    double xd1 = (v10.X - v11.X) / 2.0;
    double yd1 = (v10.Y - v11.Y) / 2.0;

    double d0 = xd0 * xd0 + yd0 * yd0;
    double d1 = xd1 * xd1 + yd1 * yd1;
    double cross = abs(xd0 * yd1 - xd1 * yd0);
    double tol = (d0 * d1) * FLT_EPSILON;

    if (cross < tol) {
        // Lines are nearly parallel
        // There are four possible minimum distance points between the lines

        double s, dx, dy, dmin = DBL_MAX;

        s = ((xd0 - xd1) * (xs0 - xs1) + (yd0 - yd1) * (ys0 - ys1)) /
            ((xd0 - xd1) * (xd0 - xd1) + (yd0 - yd1) * (yd0 - yd1));
        if (abs(s) < 1.0 - FLT_EPSILON) {
            dx = xs0 + xd0 * s - xs1 - xd1 * s;
            dy = ys0 + yd0 * s - ys1 - yd1 * s;
            dmin = fmin(dmin, dx * dx + dy * dy);

            dx = xs0 - xd0 * s - xs1 + xd1 * s;
            dy = ys0 - yd0 * s - ys1 + yd1 * s;
            dmin = fmin(dmin, dx * dx + dy * dy);
        }

        s = ((xd0 + xd1) * (xs0 - xs1) + (yd0 + yd1) * (ys0 - ys1)) /
            ((xd0 + xd1) * (xd0 + xd1) + (yd0 + yd1) * (yd0 + yd1));
        if (abs(s) < 1.0 - FLT_EPSILON) {
            dx = xs0 + xd0 * s - xs1 + xd1 * s;
            dy = ys0 + yd0 * s - ys1 + yd1 * s;
            dmin = fmin(dmin, dx * dx + dy * dy);

            dx = xs0 - xd0 * s - xs1 - xd1 * s;
            dy = ys0 - yd0 * s - ys1 - yd1 * s;
            dmin = fmin(dmin, dx * dx + dy * dy);
        }

        tol = (d0 + d1) * FLT_EPSILON;
        return (dmin < tol);
    } else { // Lines are not parallel
        double s0 = abs(xd1 * (ys0 - ys1) - yd1 * (xs0 - xs1));
        double s1 = abs(xd0 * (ys0 - ys1) - yd0 * (xs0 - xs1));
        tol = cross * (1.0 - FLT_EPSILON);

        return (s0 < tol && s1 < tol);
    }
}

bool in_triangle(Point const &p, Edge const *&e, Mesh const &m) {
    double xp = p.X;
    double yp = p.Y;

    Point const p0 = m.point(e->node());
    Point const p1 = m.point(e->next()->node());
    Point const p2 = m.point(e->prev()->node());

    double dx0p = p0.X - xp;
    double dy0p = p0.Y - yp;

    double dx1p = p1.X - xp;
    double dy1p = p1.Y - yp;

    double dx2p = p2.X - xp;
    double dy2p = p2.Y - yp;

    double dx01 = p0.X - p1.X;
    double dy01 = p0.Y - p1.Y;

    double dx12 = p1.X - p2.X;
    double dy12 = p1.Y - p2.Y;

    double dx20 = p2.X - p0.X;
    double dy20 = p2.Y - p0.Y;

    double area012 = dx01 * dy12 - dy01 * dx12;

    double tol = FLT_EPSILON * area012;

    double area01p = dx0p * dy1p - dx1p * dy0p;
    double area12p = dx1p * dy2p - dx2p * dy1p;
    double area20p = dx2p * dy0p - dx0p * dy2p;

    return (area01p > -tol && area12p > -tol && area20p > -tol);
}

void element_quality(std::vector<Edge *> &triangle, std::vector<double> &radii, std::vector<double> &quality, Mesh const &m) {
    radii.resize(0);
    quality.resize(0);

    radii.reserve(triangle.size());
    quality.reserve(triangle.size());
    for (size_t i = 0; i < triangle.size(); ++i) {
        double r = m.circumradius(triangle[i]);
        double l = m.shortest_edge_length(triangle[i]);

        radii.push_back(r);
        quality.push_back(l / r);
    }
}

void sort_permutation_ascending(std::vector<double> &values, std::vector<size_t> &index) {
    index.resize(values.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](size_t i, size_t j) { return (values[i] < values[j]); });
}

void sort_permutation_descending(std::vector<double> &values, std::vector<size_t> &index) {
    index.resize(values.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](size_t i, size_t j) { return (values[i] > values[j]); });
}