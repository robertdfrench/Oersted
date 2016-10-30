#include "Mesh.hpp"

// Edge
Edge::Edge(Curve *c, bool direction) {
    ConstraintCurve = c;
    Node = new Point(direction ? c->start() : c->end());
    Twin = this;
    Next = nullptr;
    Prev = nullptr;
    Orientation = direction;
}

Edge::Edge(Curve *c, bool direction, const Point *v) {
    ConstraintCurve = c;
    Node = v; //Should have C->point(0.0) == *(v) || c->point(1.0) == *(v)
    Twin = this;
    Next = nullptr;
    Prev = nullptr;
    Orientation = direction;
}

// Test Methods
bool Edge::is_protruding() const {
    /*
        See Chapter 1.4 of "Triangulations and Applications" by Øyvind Hjelle and Morten Dæhlen
    */
    const Point *v0 = Prev->Node;
    const Point *v1 = Node;
    const Point *v2 = Next->Node;

    double v1x = v2->X - v1->X;
    double v1y = v2->Y - v1->Y;

    double v0x = v0->X - v1->X;
    double v0y = v0->Y - v1->Y;

    double area = v1x * v0y - v1y * v0x;

    if (area > 0.0) {
        // Make sure no boundary points interior to triangle
        // Calculate barrycentric coordinates of p
        Edge * e = Next->Next;

        while (e != Prev) {
            const Point *p = e->Node;

            double v2x = v2->X - p->X;
            double v2y = v2->Y - p->Y;

            double v1x = v1->X - p->X;
            double v1y = v1->Y - p->Y;

            double v0x = v0->X - p->X;
            double v0y = v0->Y - p->Y;

            double b0 = (v0x * v1y - v0y * v1x);
            double b1 = (v1x * v2y - v1y * v2x);
            double b2 = (v2x * v0y - v2y * v0x);

            if (b0 >= 0.0 && b0 <= area && b1 >= 0.0 && b1 <= area && b2 >= 0.0 && b2 <= area) {
                return false;    // Point is interior to triangle
            } else {
                e = e->Next;
            }
        }
        return true;
    } else {
        return false;
    }
}

bool Edge::is_locally_optimal() const {
    /*
        See Chapter 3.7 of "Triangulations and Applications" by Øyvind Hjelle and Morten Dæhlen
    */

    if (ConstraintCurve != nullptr) {
        return true;
    } else {
        const Point *p3 = Node;
        const Point *p2 = Prev->Node;
        const Point *p1 = Twin->Node;
        const Point *p4 = Twin->Prev->Node;

        double v1x = p3->X - p2->X;
        double v1y = p3->Y - p2->Y;
        double v2x = p1->X - p2->X;
        double v2y = p1->Y - p2->Y;
        double v3x = p1->X - p4->X;
        double v3y = p1->Y - p4->Y;
        double v4x = p3->X - p4->X;
        double v4y = p3->Y - p4->Y;

        double d1 = sqrt(v1x * v1x + v1y * v1y);
        double d2 = sqrt(v2x * v2x + v2y * v2y);
        double d3 = sqrt(v3x * v3x + v3y * v3y);
        double d4 = sqrt(v4x * v4x + v4y * v4y);
        double tol = -(d1 * d2 * d3 * d4 * FLT_EPSILON);

        double sina = v1x * v2y - v1y * v2x;
        double sinb = v3x * v4y - v3y * v4x;
        double cosa = v1x * v2x + v1y * v2y;
        double cosb = v3x * v4x + v3y * v4y;
        double cct = sina * cosb + cosa * sinb;

        return (cct < tol ? false : true);
    }
}

bool Edge::is_valid() const {
    bool value = true;

    value = value && (this == Next->Prev);
    value = value && (this == Prev->Next);
    value = value && (this == Twin->Twin);

    return value;
}

bool Edge::is_encroached(const Point *v2) const {
    /*
    A constrained edge is encroached if a triangle and it's circumcenter lie on opposite sides of the edge.
    This is equivalent to a node being in the diameteral ball of the edge?
    This only occurs if one of the triangles attached to the edge encroaches the edge?
    This only occurs if the angle of the triangles attached to the edge has an angle opposite the edge of greater than 90 degrees.
    Using the dot product, this requires that the normalized dot product < cos(90) = 0
    */

    if (ConstraintCurve == nullptr) {
        return false;
    } else {
        const Point *v0 = base();
        const Point *v1 = tip();

        double dx0 = v0->X - v2->X;
        double dy0 = v0->Y - v2->Y;

        double dx1 = v1->X - v2->X;
        double dy1 = v1->Y - v2->Y;

        double dot = dx0 * dx1 + dy0 * dy1;
        double tol = sqrt(dx0 * dx0 + dy0 * dy0) * sqrt(dx1 * dx1 + dy1 * dy1) * FLT_EPSILON;

        return (dot < tol ? true : false);
    }
}

bool Edge::is_attached(const Point *p, Edge *&e) const {
    if (*this->tip() == *p) {
        return true;
    }

    e = const_cast<Edge *>(this);

    if (e != e->Twin) {
        e = e->Twin->Next;
        while (e != this) {
            if (*e->tip() == *p) {
                return true;
            } else if (e->Twin != e) {
                e = e->Twin->Next;
            } else {
                break;
            }
        }
    }

    if (e == e->Twin) {
        e = this->Prev;
        while (e != e->Twin) {
            if (*e->base() == *p) {
                e = e->Twin;
                return true;
            } else {
                e = e->Twin->Prev;
            }
        }
    }

    return false;
}

// Topology Altering Methods
bool Edge::swap() {
    if (ConstraintCurve == nullptr) {
        Edge * e1 = Next;
        Edge * e2 = Prev;
        Edge * e3 = Twin->Next;
        Edge * e4 = Twin->Prev;

        Node = e2->Node;
        Next = e4;
        Prev = e1;
        Mark = false;

        Twin->Node = e4->Node;
        Twin->Next = e2;
        Twin->Prev = e3;
        Twin->Mark = false;

        e1->Next = this;
        e1->Prev = e4;
        e1->Mark = false;

        e2->Next = e3;
        e2->Prev = Twin;
        e2->Mark = false;

        e3->Next = Twin;
        e3->Prev = e2;
        e3->Mark = false;

        e4->Next = e1;
        e4->Prev = this;
        e4->Mark = false;

        return true;
    } else {
        Mark = false;
        return false;
    }
}

bool Edge::recursive_swap() {
    // TODO, May need to have two different recursive swap methods, one for midpoint insertion and one for circumcenter insertion
    if (!is_locally_optimal() && swap()) {
        Edge * next = Next;
        Edge * prev = Prev;
        Edge * tnext = Twin->Next;
        Edge * tprev = Twin->Prev;

        next->recursive_swap();
        prev->recursive_swap();
        tnext->recursive_swap();
        tprev->recursive_swap();

        return true;
    } else {
        return false;
    }
}

void Edge::split_edge(std::vector<const Point *> &points, std::vector<Edge *> &edges) {
    /*
        Splits edge into two edges at the midpoint without creating any new triangles.
        Used for initial polygon refinement.
    */

    Point * p = new Point;
    points.push_back(p);

    Curve *c;
    if (ConstraintCurve != nullptr) { // Constrained Edge
        Vertex *v = new Vertex;    // #TODO: Need to manage this memory. Destructors.
        c = ConstraintCurve->split(v, 0.5);

        p->X = v->x();
        p->Y = v->y();
    } else { // Unconstrained Edge
        c = nullptr;
        *(p) = Point((Node->X + Next->Node->X) / 2.0, (Node->Y + Next->Node->Y) / 2.0);
    }

    if (this == Twin) { // Boundary Edge
        Edge * e = new Edge;
        edges.push_back(e);

        // Constraint Curve
        e->Orientation = Orientation;
        if (Orientation) {
            e->ConstraintCurve = c;
        } else {
            e->ConstraintCurve = ConstraintCurve;
            ConstraintCurve = c;
        }

        // Connectivity
        e->Node = p;
        e->Next = Next;
        e->Prev = this;
        e->Twin = e;
        e->Mark = false;

        Next->Prev = e;
        Next->Mark = false;

        Next = e;
        Mark = false;
    } else { // Interior Edge
        Edge * e0 = new Edge;
        edges.push_back(e0);

        Edge * e1 = new Edge;
        edges.push_back(e1);

        // Constraint Curve
        e0->Orientation = Orientation;
        e1->Orientation = !Orientation;
        Twin->Orientation = !Orientation;
        if (Orientation) {
            e0->ConstraintCurve = c;
        } else {
            e0->ConstraintCurve = ConstraintCurve;
            ConstraintCurve = c;
        }
        Twin->ConstraintCurve = e0->ConstraintCurve;
        e1->ConstraintCurve = ConstraintCurve;

        // Connectivity
        e0->Node = p;
        e0->Next = Next;
        e0->Prev = this;
        e0->Twin = Twin;
        e0->Mark = false;

        e1->Node = p;
        e1->Next = Twin->Next;
        e1->Prev = Twin;
        e1->Twin = this;
        e1->Mark = false;

        if (Twin->Next != nullptr) {
            Twin->Next->Prev = e1;
            Twin->Next->Mark = false;
        }

        Twin->Next = e1;
        Twin->Twin = e0;
        Twin->Mark = false;

        if (Next != nullptr) {
            Next->Prev = e0;
            Next->Mark = false;
        }

        Next = e0;
        Twin = e1;
        Mark = false;
    }
}

// Calculation Methods
Point Edge::circumcenter() const {
    double xa = Node->X;
    double ya = Node->Y;
    double xb = Prev->Node->X - xa;
    double yb = Prev->Node->Y - ya;
    double xc = Next->Node->X - xa;
    double yc = Next->Node->Y - ya;

    double d = 2.0 * (xb * yc - yb * xc);
    double db = xb * xb + yb * yb;
    double dc = xc * xc + yc * yc;
    double x = (yc * db - yb * dc) / d + xa;
    double y = (xb * dc - xc * db) / d + ya;

    return Point{x, y};
}

double Edge::circumradius() const {
    double xa = Node->X;
    double ya = Node->Y;
    double xb = Prev->Node->X - xa;
    double yb = Prev->Node->Y - ya;
    double xc = Next->Node->X - xa;
    double yc = Next->Node->Y - ya;
    xa = xb - xc;
    ya = yb - yc;

    double den = 2.0 * abs(xb * yc - yb * xc);
    double num = sqrt(xa * xa + ya * ya) * sqrt(xb * xb + yb * yb) * sqrt(xc * xc + yc * yc);

    return num / den;
}

double Edge::length() const {
    double dx = Node->X - Next->Node->X;
    double dy = Node->Y - Next->Node->Y;

    return sqrt(dx * dx + dy * dy);
}

double Edge::shortest_edge_length() const {
    double x0 = Node->X;
    double y0 = Node->Y;
    double x1 = Next->Node->X;
    double y1 = Next->Node->Y;
    double x2 = Prev->Node->X;
    double y2 = Prev->Node->Y;

    double dx = x0 - x1;
    double dy = y0 - y1;
    double dl = dx * dx + dy * dy;

    dx = x1 - x2;
    dy = y1 - y2;
    dl = fmin(dl, dx * dx + dy * dy);

    dx = x2 - x0;
    dy = y2 - y0;
    dl = fmin(dl, dx * dx + dy * dy);

    return sqrt(dl);
}

// Auxillary Method
void Edge::recursive_mark() {
    if (Mark && Next->Mark && Prev->Mark) {
        Next->Mark = false;
        Prev->Mark = false;
    }
}

// Friend Utility
bool are_intersecting(const Edge *e0, const Edge *e1) {
    // #TODO, Make more detailed return type enumeration
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

// External Utility
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