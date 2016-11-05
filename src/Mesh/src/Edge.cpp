#include "Mesh.hpp"

Edge::Edge(Curve *c, bool direction) {
    ConstraintCurve = c;
    Node = new Point(direction ? c->start() : c->end());
    Twin = this;
    Next = nullptr;
    Prev = nullptr;
    Orientation = direction;
}

Edge::Edge(Curve *c, bool direction, Point const *v) {
    ConstraintCurve = c;
    Node = v; //Should have C->point(0.0) == *(v) || c->point(1.0) == *(v)
    Twin = this;
    Next = nullptr;
    Prev = nullptr;
    Orientation = direction;
}

bool Edge::is_protruding() const {
    /*
        See Chapter 1.4 of "Triangulations and Applications" by Øyvind Hjelle and Morten Dæhlen
    */
    Point const *v0 = Prev->Node;
    Point const *v1 = Node;
    Point const *v2 = Next->Node;

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
            Point const *p = e->Node;

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
        Point const *p3 = Node;
        Point const *p2 = Prev->Node;
        Point const *p1 = Twin->Node;
        Point const *p4 = Twin->Prev->Node;

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

bool Edge::is_encroached(Point const *v2) const {
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
        Point const *v0 = base();
        Point const *v1 = tip();

        double dx0 = v0->X - v2->X;
        double dy0 = v0->Y - v2->Y;

        double dx1 = v1->X - v2->X;
        double dy1 = v1->Y - v2->Y;

        double dot = dx0 * dx1 + dy0 * dy1;
        double tol = sqrt(dx0 * dx0 + dy0 * dy0) * sqrt(dx1 * dx1 + dy1 * dy1) * FLT_EPSILON;

        return (dot < tol ? true : false);
    }
}

bool Edge::is_attached(Point const *p, Edge *&e) const {
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
        Edge * enext = Next;
        Edge * eprev = Prev;
        Edge * tnext = Twin->Next;
        Edge * tprev = Twin->Prev;

        enext->recursive_swap();
        eprev->recursive_swap();
        tnext->recursive_swap();
        tprev->recursive_swap();

        return true;
    } else {
        return false;
    }
}

void Edge::split_edge(std::vector<Point const *> &points, std::vector<Edge *> &edges) {
    /*
        Splits edge into two edges at the midpoint without creating any new triangles.
        Used for initial polygon refinement.
    */

    // TODO: Refactor into split_constrainted_edge

    Point * p = new Point;
    points.push_back(p);

    Curve *c;
    if (ConstraintCurve != nullptr) { // Constrained Edge
        Vertex *v = new Vertex;    // TODO: Need to manage this memory. Destructors.
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
    } else {
        throw std::exception(); // Function should only be called on
    }
}

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

void Edge::recursive_mark() {
    if (Mark && Next->Mark && Prev->Mark) {
        Next->Mark = false;
        Prev->Mark = false;
    }
}