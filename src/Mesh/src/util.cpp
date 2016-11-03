#include "Mesh.hpp"

bool are_intersecting(const Edge *e0, const Edge *e1) {
    // TODO, Make more detailed return type enumeration
    if (e0->ConstraintCurve != nullptr && e0->ConstraintCurve == e1->ConstraintCurve) {
        return false;
    }

    const Point *v00 = e0->base();
    const Point *v01 = e0->tip();
    const Point *v10 = e1->base();
    const Point *v11 = e1->tip();

    double xs0 = (v00->X + v01->X) / 2.0;
    double ys0 = (v00->Y + v01->Y) / 2.0;
    double xs1 = (v10->X + v11->X) / 2.0;
    double ys1 = (v10->Y + v11->Y) / 2.0;

    double xd0 = (v00->X - v01->X) / 2.0;
    double yd0 = (v00->Y - v01->Y) / 2.0;
    double xd1 = (v10->X - v11->X) / 2.0;
    double yd1 = (v10->Y - v11->Y) / 2.0;

    double d0 = xd0 * xd0 + yd0 * yd0;
    double d1 = xd1 * xd1 + yd1 * yd1;
    double cross = abs(xd0 * yd1 - xd1 * yd0);
    double tol = (d0 * d1) * FLT_EPSILON;

    if (cross < tol) {
        /*
            Lines are nearly parallel
            There are four possible minimum distance points between the lines
        */
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
    } else {
        // Lines are not parallel
        double s0 = abs(xd1 * (ys0 - ys1) - yd1 * (xs0 - xs1));
        double s1 = abs(xd0 * (ys0 - ys1) - yd0 * (xs0 - xs1));
        tol = cross * (1.0 - FLT_EPSILON);

        return (s0 < tol && s1 < tol);
    }
}

void circumradius(std::vector<Edge *> &triangles, std::vector<double> &radii) {
    radii.resize(0);
    radii.reserve(triangles.size());
    for (size_t i = 0; i < triangles.size(); ++i) {
        radii.push_back(triangles[i]->circumradius());
    }
}

void shortest_edge_length(std::vector<Edge *> &triangle, std::vector<double> &length) {
    length.resize(0);
    length.reserve(triangle.size());
    for (size_t i = 0; i < triangle.size(); ++i) {
        length.push_back(triangle[i]->shortest_edge_length());
    }
}

void element_quality(std::vector<Edge *> &triangle, std::vector<double> &radii, std::vector<double> &quality) {
    radii.resize(0);
    quality.resize(0);

    radii.reserve(triangle.size());
    quality.reserve(triangle.size());
    for (size_t i = 0; i < triangle.size(); ++i) {
        double r = triangle[i]->circumradius();
        double l = triangle[i]->shortest_edge_length();

        radii.push_back(r);
        quality.push_back(l / r);
    }
}

void sort_permutation(std::vector<double> &values, std::vector<size_t> &index) {
    index.resize(values.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](size_t i, size_t j) { return (values[i] < values[j]); });
}