#include "Mesh.hpp"

Mesh::Mesh(Sketch &sketch) {
    Boundary = sketch.boundary();

    for (size_t i = 0; i != sketch.size_curves(); ++i) {
        auto c = sketch.curve(i);
        if(!(c->ForConstruction)) {
            Curves.push_back(c);
        }
    }

    for (size_t i = 0; i != sketch.size_contours(); ++i) {
        Contours.push_back(sketch.contour(i));
    }
}

bool Mesh::are_intersecting(Edge const *e0, Edge const *e1) const {
    // TODO, Make more detailed return type enumeration
    if (e0->ConstraintCurve != nullptr && e0->ConstraintCurve == e1->ConstraintCurve) {
        return false;
    }

    Point const v00 = Points[e0->base()];
    Point const v01 = Points[e0->tip(*this)];
    Point const v10 = Points[e1->base()];
    Point const v11 = Points[e1->tip(*this)];

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

bool Mesh::edges_are_valid() const {
    bool result = true;

    for (auto e : Edges) {
        if (e->Self != Edges[e->Next]->Prev) {
            result = false;
            break;
        }
        if (e->Self != Edges[e->Prev]->Next) {
            result = false;
            break;
        }
        if (e->Self != Edges[e->Twin]->Twin) {
            result = false;
            break;
        }

        if ((e->Twin != e->Self)) {
            if (e->node() != Edges[Edges[e->Twin]->Next]->node()) {
                result = false;
                break;
            }
            if (e->constraint_curve() != Edges[e->twin()]->constraint_curve()) {
                result = false;
                break;
            }
            if (e->is_constrained()) {
                if (e->orientation() == Edges[e->twin()]->orientation()) {
                    result = false;
                    break;
                }
            }

            if (e->node() == Edges[e->twin()]->node()) {
                result = false;
                break;
            }
        }

        if (e->is_constrained()) {
            if (e->orientation()) {
                if (Points[e->base()] != *e->constraint_curve()->start()) {
                    result = false;
                    break;
                }
                if (Points[e->tip(*this)] != *e->constraint_curve()->end()) {
                    result = false;
                    break;
                }
            } else {
                if (Points[e->base()] != *e->constraint_curve()->end()) {
                    result = false;
                    break;
                }
                if (Points[e->tip(*this)] != *e->constraint_curve()->start()) {
                    result = false;
                    break;
                }
            }
        }
    }

    return result;
}

bool Mesh::find_attached(Edge *&e_out, Point const p) const {
    if (Points[e_out->tip(*this)] == p) {
        return true;
    }

    Edge * e_in = e_out;

    if (e_out->Self != e_out->Twin) {
        e_out = Edges[Edges[e_out->Twin]->Next];
        while (e_out->Self != e_in->Self) {
            if (Points[e_out->tip(*this)] == p) {
                return true;
            } else if (e_out->Twin != e_out->Self) {
                e_out = Edges[Edges[e_out->Twin]->Next];
            } else {
                break;
            }
        }
    }

    if (e_out->Self == e_out->Twin) {
        e_out = Edges[e_in->Prev];
        while (e_out->Self != e_out->Twin) {
            if (Points[e_out->base()] == p) {
                e_out = Edges[e_out->Twin];
                return true;
            } else {
                e_out = Edges[Edges[e_out->Twin]->Prev];
            }
        }
    }

    return false;
}

bool Mesh::refine() {
    std::vector<double> radii;
    std::vector<double> quality;
    std::vector<size_t> index;

    // #TODO, Loop until quality is satisfied
    element_quality(Triangles, radii, quality);
    sort_permutation_ascending(quality, index);
    size_t N = Triangles.size();

    refine_once(index, radii, quality);
    size_t M = Triangles.size();

    while (M > N) {
        N = M;
        element_quality(Triangles, radii, quality);
        sort_permutation_ascending(quality, index);
        refine_once(index, radii, quality);
        M = Triangles.size();
    }

    return edges_are_valid();
}

bool Mesh::refine_once() {
    std::vector<double> radii;
    std::vector<double> quality;
    std::vector<size_t> index;

    element_quality(Triangles, radii, quality);
    sort_permutation_ascending(quality, index);
    //sort_permutation_descending(radii, index);
    refine_once(index, radii, quality);

    return edges_are_valid();
}

bool Mesh::in_triangle(Point const p, Edge const *e) const {
    double xp = p.X;
    double yp = p.Y;

    Point const p0 = Points[e->Node];
    Point const p1 = Points[Edges[e->Next]->Node];
    Point const p2 = Points[Edges[e->Prev]->Node];

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

bool Mesh::is_encroached(Edge const *e, Point const p) const {
    /*
        A constrained edge is encroached if a triangle and it's circumcenter lie on opposite sides of the edge.
        This is equivalent to a node being in the diameteral ball of the edge?
        This only occurs if one of the triangles attached to the edge encroaches the edge?
        This only occurs if the angle of the triangles attached to the edge has an angle opposite the edge of greater than 90 degrees.
        Using the dot product, this requires that the normalized dot product < cos(90) = 0
    */

    if (!(e->is_constrained())) {
        return false;
    } else {
        Point const p0 = Points[e->base()];
        Point const p1 = Points[e->tip(*this)];

        double dx0 = p0.X - p.X;
        double dy0 = p0.Y - p.Y;

        double dx1 = p1.X - p.X;
        double dy1 = p1.Y - p.Y;

        double dot = dx0 * dx1 + dy0 * dy1;
        double tol = sqrt(dx0 * dx0 + dy0 * dy0) * sqrt(dx1 * dx1 + dy1 * dy1) * FLT_EPSILON;

        return (dot < tol);
    }
}

bool Mesh::is_locally_optimal(Edge const *e) const {
    /*
        See Chapter 3.7 of "Triangulations and Applications" by Øyvind Hjelle and Morten Dæhlen
    */

    if (e->is_constrained()) {
        return true;
    } else {
        Point const p3 = Points[e->Node];
        Point const p2 = Points[Edges[e->Prev]->Node];
        Point const p1 = Points[Edges[e->Twin]->Node];
        Point const p4 = Points[Edges[Edges[e->Twin]->Prev]->Node];

        double v1x = p3.X - p2.X;
        double v1y = p3.Y - p2.Y;
        double v2x = p1.X - p2.X;
        double v2y = p1.Y - p2.Y;
        double v3x = p1.X - p4.X;
        double v3y = p1.Y - p4.Y;
        double v4x = p3.X - p4.X;
        double v4y = p3.Y - p4.Y;

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

        return cct > tol;
    }
}

bool Mesh::is_protruding(Edge const *e) const {
    //
    //    See Chapter 1.4 of "Triangulations and Applications" by Øyvind Hjelle and Morten Dæhlen
    //

    Point const p0 = Points[Edges[e->Prev]->Node];
    Point const p1 = Points[e->Node];
    Point const p2 = Points[Edges[e->Next]->Node];

    double v1x = p2.X - p1.X;
    double v1y = p2.Y - p1.Y;

    double v0x = p0.X - p1.X;
    double v0y = p0.Y - p1.Y;

    double area = v1x * v0y - v1y * v0x;

    if (area > 0.0) {
        // Make sure no boundary points interior to triangle
        // Calculate barrycentric coordinates of p
        Edge * next = Edges[Edges[e->Next]->Next];

        while (next->Self != e->Prev) {
            Point const p4 = Points[next->Node];

            double v2x = p2.X - p4.X;
            double v2y = p2.Y - p4.Y;

            v1x = p1.X - p4.X;
            v1y = p1.Y - p4.Y;

            v0x = p0.X - p4.X;
            v0y = p0.Y - p4.Y;

            double b0 = (v0x * v1y - v0y * v1x);
            double b1 = (v1x * v2y - v1y * v2x);
            double b2 = (v2x * v0y - v2y * v0x);

            if (b0 >= 0.0 && b0 <= area && b1 >= 0.0 && b1 <= area && b2 >= 0.0 && b2 <= area) {
                return false;    // Point is interior to triangle
            } else {
                next = Edges[next->Next];
            }
        }
        return true;
    } else {
        return false;
    }
}

bool Mesh::is_valid(Edge const *e) const {
    bool value = true;

    value = value && (e->Self == Edges[e->Next]->Prev);
    value = value && (e->Self == Edges[e->Prev]->Next);
    value = value && (e->Self == Edges[e->Twin]->Twin);

    return value;
}

bool Mesh::recursive_swap(Edge *e) const {
    // TODO, May need to have two different recursive swap methods, one for midpoint insertion and one for circumcenter insertion
    if (!is_locally_optimal(e) && swap(e)) {
        Edge * enext = Edges[e->Next];
        Edge * eprev = Edges[e->Prev];
        Edge * tnext = Edges[Edges[e->Twin]->Next];
        Edge * tprev = Edges[Edges[e->Twin]->Prev];

        recursive_swap(enext);
        recursive_swap(eprev);
        recursive_swap(tnext);
        recursive_swap(tprev);

        return true;
    } else {
        return false;
    }
}

bool Mesh::swap(Edge *&e0) const {
    if (!e0->is_constrained()) {
        Edge *e1 = Edges[e0->Next];
        Edge *e2 = Edges[e0->Prev];
        Edge *e3 = Edges[Edges[e0->Twin]->Next];
        Edge *e4 = Edges[Edges[e0->Twin]->Prev];
        Edge *twin = Edges[e0->Twin];

        e0->Node = e2->Node;
        e0->Next = e4->Self;
        e0->Prev = e1->Self;
        e0->Mark = false;

        twin->Node = e4->Node;
        twin->Next = e2->Self;
        twin->Prev = e3->Self;
        twin->Mark = false;

        e1->Next = e0->Self;
        e1->Prev = e4->Self;
        e1->Mark = false;

        e2->Next = e3->Self;
        e2->Prev = e0->Twin;
        e2->Mark = false;

        e3->Next = e0->Twin;
        e3->Prev = e2->Self;
        e3->Mark = false;

        e4->Next = e1->Self;
        e4->Prev = e0->Self;
        e4->Mark = false;

        return true;
    } else {
        e0->Mark = false;
        return false;
    }
}

double Mesh::circumradius(Edge const *e) const {
    Point const node = Points[e->Node];
    Point const prev = Points[Edges[e->Prev]->Node];
    Point const next = Points[Edges[e->Next]->Node];

    double xa = node.X;
    double ya = node.Y;
    double xb = prev.X - xa;
    double yb = prev.Y - ya;
    double xc = next.X - xa;
    double yc = next.Y - ya;
    xa = xb - xc;
    ya = yb - yc;

    double den = 2.0 * abs(xb * yc - yb * xc);
    double num = sqrt(xa * xa + ya * ya) * sqrt(xb * xb + yb * yb) * sqrt(xc * xc + yc * yc);

    return num / den;
}

double Mesh::length(Edge const *e) const {
    Point const base = Points[e->base()];
    Point const tip = Points[e->tip(*this)];

    double dx = base.X - tip.X;
    double dy = base.Y - tip.Y;

    return sqrt(dx * dx + dy * dy);
}

double Mesh::shortest_edge_length(Edge const *e) const {
    Point const node = Points[e->Node];
    Point const prev = Points[Edges[e->Prev]->Node];
    Point const next = Points[Edges[e->Next]->Node];

    double x0 = node.X;
    double y0 = node.Y;
    double x1 = next.X;
    double y1 = next.Y;
    double x2 = prev.X;
    double y2 = prev.Y;

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

size_t Mesh::num_edges() const {
    size_t count = 0;
    for (size_t i = 0; i < Edges.size(); ++i) {
        count += (Edges[i]->Self == Edges[i]->Twin ? 2 : 1);
    }
    count /= 2;

    return count;
}

void Mesh::create() {
    create_boundary_polygon();
    triangulate_boundary_polygon();
    get_triangles();

    insert_internal_boundaries();
    get_triangles();
}

void Mesh::create_boundary_polygon() {
    // Create input edges
    Edges.reserve(Boundary->size());
    Points.reserve(Boundary->size());
    for (size_t i = 0; i != Boundary->size(); ++i) {
        Curve * c = Boundary->curve(i)->clone(); //clone() to prevent alteration of input Contour when Edge is split
        bool dir = Boundary->orientation(i);
        Edges.push_back(new Edge(c, dir, Points.size(), Edges.size())); // TODO: Write new_edge function
        Edges.back()->Twin = Edges.back()->Self;
        Points.push_back(Point(dir ? c->start() : c->end()));
    }

    // Set Next/Prev edges from Boundary
    size_t Ne = 0;
    for (size_t i = 0; i != Edges.size(); ++i) {
        size_t j = (i + 1) % Edges.size();
        Edges[i]->Next = Edges[j]->Self;
        Edges[j]->Prev = Edges[i]->Self;
    }

    // TODO: Length Based Refinement of Boundary curves

    // Some edges may be intersect due to discretization error
    // If two edges intersect, split the longest edge
    bool any_split = true;
    while (any_split) {
        any_split = false;
        for (size_t i = 0; i != Edges.size(); ++i) {
            for (size_t j = i + 1; j != Edges.size(); ++j) {
                if (are_intersecting(Edges[i], Edges[j])) {
                    any_split = true;
                    if (length(Edges[i]) > length(Edges[j])) {
                        split_edge(Edges[i]);
                    } else {
                        split_edge(Edges[j]);
                    }
                }
            }
        }
    }
}

void Mesh::delete_me() {
    for (auto &ed : Edges) {
        delete ed;
    }
    Edges.clear();

    Points.clear();
}

void Mesh::element_quality(std::vector<Edge *> &triangle, std::vector<double> &radii, std::vector<double> &quality) {
    radii.resize(0);
    quality.resize(0);

    radii.reserve(triangle.size());
    quality.reserve(triangle.size());
    for (size_t i = 0; i < triangle.size(); ++i) {
        double r = circumradius(triangle[i]);
        double l = shortest_edge_length(triangle[i]);

        radii.push_back(r);
        quality.push_back(l / r);
    }
}

void Mesh::get_triangles() {
    Triangles.resize(0);
    Triangles.reserve(2 * num_points());

    mark_triangles();
    for (size_t i = 0; i < Edges.size(); ++i) {
        if (Edges[i]->Mark) {
            Triangles.push_back(Edges[i]);
        }
    }
}

void Mesh::insert_internal_boundaries() {
    /*
        For each constraint curve :
            Add constraint curve to queue
            Insert endpoints if they do not exist.
                Will have to keep attempting insertion until no existing edge is encroached

            While the queue is not empty :
                Orbit the start point of the last curve in the queue :
                    If the end point is attached to the start point by some edge :
                        Set Edge->ConstraintCurve and->Orientation properties
                        Pop last curve from queue
                    Else
                        Split the last curve and add a new curve to the end of the queue

            Repeat until no edge is split :
                Orbit each edge, checking for encroachment by attached verticies
                If encroached, split edge
    */

    // Find interior curves
    std::vector<Curve *> interior;
    for (auto c : Curves) {
        bool on_exterior = false;
        for (size_t i = 0; i != Boundary->size(); ++i) {
            if (c == Boundary->curve(i)) {
                on_exterior = true;
                break;
            }
        }

        if (!on_exterior) {
            interior.push_back(c->clone());
        }
    }

    // Insert interior curve end points
    for (size_t i = 0; i != interior.size(); ++i) {
        // Insert start point
        Point p = interior[i]->start();
        LocateTriangleResult result = locate_triangle(p);
        if (result == LocateTriangleResult::Interior) {
            while (insert_point(p) == InsertPointResult::Midpoint);
        }

        // Insert end point
        p = interior[i]->end();
        result = locate_triangle(p);
        if (result == LocateTriangleResult::Interior) {
            while (insert_point(p) == InsertPointResult::Midpoint);
        }
    }

    // Insert interior curve midpoints until constraints are naturally satisfied
    std::vector<Curve *> queue;
    for (size_t i = 0; i != interior.size(); ++i) {
        Edge * e = Edges.back();
        queue.push_back(interior[i]);
        while (queue.size() != 0) {
            Point p0 = queue.back()->start();
            Point p1 = queue.back()->end();

            LocateTriangleResult result = locate_triangle(p0, e);

            if (result != LocateTriangleResult::Point) {
                throw std::exception();
            }

            if (find_attached(e, p1)) {
                e->ConstraintCurve = queue.back();
                e->Orientation = true;

                Edges[e->Twin]->ConstraintCurve = queue.back();
                Edges[e->Twin]->Orientation = false;

                queue.pop_back();
            } else {
                Vertex *v = new Vertex; // TODO: Manage memory
                queue.push_back(queue.back()->split(v, 0.5));

                Point const p = queue.back()->start();

                while (insert_point(p) == InsertPointResult::Midpoint);
            }
        }
    }

    split_encroached_edges();
}

void Mesh::mark_triangles() {
    for (size_t i = 0; i < Edges.size(); ++i) {
        Edges[i]->Mark = true;
    }

    for (size_t i = 0; i < Edges.size(); ++i) {
        recursive_mark(Edges[i]);
    }
}

void Mesh::recursive_mark(Edge *e) const {
    if (e->Mark && Edges[e->Next]->Mark && Edges[e->Prev]->Mark) {
        Edges[e->Next]->Mark = false;
        Edges[e->Prev]->Mark = false;
    }
}

void Mesh::refine_once(std::vector<size_t> index, std::vector<double> radii, std::vector<double> quality) {
    for (size_t i = 0; i < Triangles.size(); ++i) {
        size_t j = index[i];
        if ((Triangles[j]->Mark) && ((radii[j] > MaximumElementSize) || (radii[j] > MinimumElementSize && quality[j] < MinimumElementQuality))) {
            insert_circumcenter(Triangles[j]);
        }
    }
    get_triangles();
}

void Mesh::save_as(std::string path, std::string file_name) const {
    /*
        This is a stub for visualization
    */

    if (!boost::filesystem::exists(path)) {
        boost::filesystem::create_directories(path);
    }

    std::fstream fs;
    fs.open(path + file_name + ".oeme", std::fstream::out);

    for (size_t i = 0; i < Edges.size(); i++) {
        Point const v0 = Points[Edges[i]->Node];
        Point const v1 = Points[Edges[Edges[i]->Next]->Node];
        Point const v2 = Points[Edges[Edges[Edges[i]->Next]->Next]->Node];
        fs << v0.X << ',' << v1.X << ',' << v2.X << ',' << v0.Y << ',' << v1.Y << ',' << v2.Y << '\n';
    }

    fs.close();
}

void Mesh::sort_permutation_ascending(std::vector<double> &values, std::vector<size_t> &index) const {
    index.resize(values.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](size_t i, size_t j) { return (values[i] < values[j]); });
}

void Mesh::sort_permutation_descending(std::vector<double> &values, std::vector<size_t> &index) const {
    index.resize(values.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](size_t i, size_t j) { return (values[i] > values[j]); });
}

void Mesh::split_edge(Edge *e) {
    /*
        Splits edge into two edges at the midpoint without creating any new triangles.
        Used for initial polygon refinement.
    */

    // TODO: Refactor into split_constrained_edge

    Curve *c;
    if (e->is_constrained()) {
        Vertex *v = new Vertex;    // TODO: Need to manage this memory.
        c = e->ConstraintCurve->split(v, 0.5);

        Points.push_back(Point(v->x(),v->y()));
    } else { // Unconstrained Edge
        c = nullptr;
        Point const base = Points[e->base()];
        Point const tip = Points[e->tip(*this)]; // TODO:

        Points.push_back(Point((base.X + tip.X) / 2.0, (base.Y + tip.Y) / 2.0));
    }

    if (e->Self == e->Twin) { // Boundary Edge, TODO: write bool Edge::is_boundary()
        Edge * newe = new Edge;
        add_edge(newe);

        // Constraint Curve
        newe->Orientation = e->Orientation;
        if (e->Orientation) {
            newe->ConstraintCurve = c;
        } else {
            newe->ConstraintCurve = e->ConstraintCurve;
            e->ConstraintCurve = c;
        }

        // Connectivity
        newe->Node = Points.size() - 1;
        newe->Next = e->Next;
        newe->Prev = e->Self;
        newe->Twin = newe->Self;
        newe->Mark = false;

        Edges[e->Next]->Prev = newe->Self;
        Edges[e->Next]->Mark = false;

        e->Next = newe->Self;
        e->Mark = false;
    } else {
        throw std::exception(); // Function should only be called on constrained edges
    }
}

void Mesh::split_encroached_edges() {
    bool any_split = true;
    while (any_split) {
        any_split = false;
        for (size_t i = 0; i != Edges.size(); ++i) {
            if (Edges[i]->is_constrained()) {
                if (is_encroached(Edges[i], Points[Edges[Edges[i]->Prev]->base()])) {
                    any_split = true;
                    insert_midpoint(Edges[i]);
                }
            }
        }
    }
}

void Mesh::triangulate_boundary_polygon() {
    Edges.reserve(3 * num_points());

    Edge * ep = Edges[0];
    while (ep->Self != Edges[Edges[ep->Next]->Next]->Next) {
        if (is_protruding(ep)) {
            Edge * e0 = new Edge;
            add_edge(e0);

            Edge * e1 = new Edge;
            add_edge(e1);

            // Edge of new triangle
            e0->Node = Edges[ep->Next]->Node;
            e0->Prev = ep->Self;
            e0->Next = ep->Prev;
            e0->Twin = e1->Self;

            // Twin edge, part of new polygonal boundary
            e1->Node = Edges[ep->Prev]->Node;
            e1->Next = ep->Next;
            e1->Prev = Edges[ep->Prev]->Prev;
            e1->Twin = e0->Self;

            // Update polygonal boundary
            Edges[ep->Next]->Prev = e1->Self;
            ep->Next = e0->Self;
            Edges[Edges[ep->Prev]->Prev]->Next = e1->Self;
            Edges[ep->Prev]->Prev = e0->Self;

            // Next edge
            ep = Edges[e1->Next];
        } else {
            ep = Edges[ep->Next];
        }
    }

    // Edge swap to make triangulation Delaunay
    for (size_t i = 0; i < Edges.size(); ++i) {
        recursive_swap(Edges[i]);
    }

    split_encroached_edges();
}

LocateTriangleResult Mesh::locate_triangle(Point const p, Edge *&e) const {
    Edge const *ec = e;

    LocateTriangleResult result = locate_triangle(p, ec);

    e = const_cast<Edge *>(ec);

    return result;
}

LocateTriangleResult Mesh::locate_triangle(Point const p, Edge const *&e) const {
    double xp = p.X;
    double yp = p.Y;

    Point p0 = Points[e->Node];
    Point p1 = Points[Edges[e->Next]->Node];
    Point p2 = Points[Edges[e->Prev]->Node];

    if (p == p0) {
        return LocateTriangleResult::Point;
    } else if (p == p1) {
        e = Edges[e->Next];
        return LocateTriangleResult::Point;
    } else if (p == p2) {
        e = Edges[e->Prev];
        return LocateTriangleResult::Point;
    }

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

    double tol = FLT_EPSILON * (dx20 * dy01 - dy20 * dx01);

    double area01 = dx0p * dy1p - dx1p * dy0p;
    double area12 = dx1p * dy2p - dx2p * dy1p;
    double area20 = dx2p * dy0p - dx0p * dy2p;

    if (area01 > tol && area12 > tol && area20 > tol) {
        return LocateTriangleResult::Interior;
    } else if (area01 < -tol && e->Self != e->Twin) {
        e = Edges[e->Twin];
        p2 = p1;
        p1 = p0;
        p0 = p2;

        dx2p = dx1p;
        dx1p = dx0p;
        dx0p = dx2p;

        dy2p = dy1p;
        dy1p = dy0p;
        dy0p = dy2p;

        dx01 = -dx01;
        dy01 = -dy01;
    } else if (area12 < -tol && e->Next != Edges[e->Next]->Twin) {
        e = Edges[Edges[e->Next]->Twin];
        p0 = p2;

        dx0p = dx2p;
        dy0p = dy2p;

        dx01 = -dx12;
        dy01 = -dy12;
    } else if (area20 < -tol && e->Prev != Edges[e->Prev]->Twin) {
        e = Edges[Edges[e->Prev]->Twin];
        p1 = p2;

        dx1p = dx2p;
        dy1p = dy2p;

        dx01 = -dx20;
        dy01 = -dy20;
    } else if (area01 > -tol && area12 > tol && area20 > tol) {
        e = Edges[e->Twin];
        return LocateTriangleResult::Interior;
    } else if (area01 > tol && area12 > -tol && area20 > tol) {
        e = Edges[Edges[e->Next]->Twin];
        return LocateTriangleResult::Interior;
    } else if (area01 > tol && area12 > tol && area20 > -tol) {
        e = Edges[Edges[e->Prev]->Twin];
        return LocateTriangleResult::Interior;
    } else if (area01 < -tol) {
        return LocateTriangleResult::Exterior;
    } else if (area12 < -tol) {
        e = Edges[e->Next];
        return LocateTriangleResult::Exterior;
    } else if (area20 < -tol) {
        e = Edges[e->Prev];
        return LocateTriangleResult::Exterior;
    } else {
        throw std::exception();
    }

    while (true) { // After first iteration, area01 > 0
        p2 = Points[Edges[e->Prev]->Node];

        if (p == p2) {
            e = Edges[e->Prev];
            return LocateTriangleResult::Point;
        }

        dx2p = p2.X - xp;
        dy2p = p2.Y - yp;

        dx12 = p1.X - p2.X;
        dy12 = p1.Y - p2.Y;

        dx20 = p2.X - p0.X;
        dy20 = p2.Y - p0.Y;

        tol = FLT_EPSILON * (dx20 * dy01 - dy20 * dx01);

        area12 = dx1p * dy2p - dx2p * dy1p;
        area20 = dx2p * dy0p - dx0p * dy2p;

        if (area12 > tol && area20 > tol) {
            return LocateTriangleResult::Interior;
        } else if (area12 < -tol && e->Next != Edges[e->Next]->Twin) {
            e = Edges[Edges[e->Next]->Twin];
            p0 = p2;

            dx0p = dx2p;
            dy0p = dy2p;

            dx01 = -dx12;
            dy01 = -dy12;
            continue;
        } else if (area20 < -tol && e->Prev != Edges[e->Prev]->Twin) {
            e = Edges[Edges[e->Prev]->Twin];
            p1 = p2;

            dx1p = dx2p;
            dy1p = dy2p;

            dx01 = -dx20;
            dy01 = -dy20;
            continue;
        } else if (area12 > -tol && area20 > tol) {
            e = Edges[Edges[e->Next]->Twin];
            return LocateTriangleResult::Interior;
        } else if (area12 > tol && area20 > -tol) {
            e = Edges[Edges[e->Prev]->Twin];
            return LocateTriangleResult::Interior;
        } else if (area12 < -tol) {
            e = Edges[e->Next];
            return LocateTriangleResult::Exterior;
        } else if (area20 < -tol) {
            e = Edges[e->Prev];
            return LocateTriangleResult::Exterior;
        } else {
            throw std::exception();
        }
    }
}

InsertPointResult Mesh::insert_circumcenter(Edge *t) {
    return insert_point(circumcenter(t), t);
}

InsertPointResult Mesh::insert_midpoint(Edge *e) {
    /*
        Splits edge into two edges and creates two new triangles.
    */

    Curve *c;
    if (e->is_constrained()) { // Constrained Edge
        Vertex *v = new Vertex; // TODO: Manage memory
        c = e->ConstraintCurve->split(v, 0.5);
        Points.push_back(Point(v->x(), v->y()));
    } else { // Unconstrained Edge
        c = nullptr;
        Point const base = Points[e->base()];
        Point const tip = Points[e->tip(*this)];
        Points.push_back(Point((base.X + tip.X) / 2.0, (base.Y + tip.Y) / 2.0));
    }

    if (e->Self == e->Twin) { // Boundary Edge
        Edge * e0 = new Edge;
        add_edge(e0);

        Edge * e1 = new Edge;
        add_edge(e1);

        Edge * e2 = new Edge;
        add_edge(e2);

        // Save edges for swapping
        Edge * next = Edges[e->Next];
        Edge * prev = Edges[e->Prev];

        // Handle constraint curves
        e2->Orientation = e->Orientation;
        if (e->Orientation) {
            e2->ConstraintCurve = c;
        } else {
            e2->ConstraintCurve = e->ConstraintCurve;
            e->ConstraintCurve = c;
        }

        // Construct edges
        e0->Node = Edges[e->Prev]->Node;
        e0->Next = e2->Self;
        e0->Prev = e->Next;
        e0->Twin = e1->Self;
        e0->Mark = false;

        e1->Node = Points.size() - 1; // TODO: write add_point() method
        e1->Next = e->Prev;
        e1->Prev = e->Self;
        e1->Twin = e0->Self;
        e1->Mark = false;

        e2->Node = Points.size() - 1;
        e2->Next = e->Next;
        e2->Prev = e0->Self;
        e2->Twin = e2->Self;
        e2->Mark = false;

        Edges[e->Next]->Next = e0->Self;
        Edges[e->Next]->Prev = e2->Self;
        Edges[e->Next]->Mark = false;

        Edges[e->Prev]->Prev = e1->Self;
        Edges[e->Prev]->Mark = false;

        e->Next = e1->Self;
        e->Mark = false;

        // Recursive swap
        recursive_swap(e0);
        recursive_swap(next);
        recursive_swap(prev);
    } else { // Interior Edge
        Edge * e0 = new Edge;
        add_edge(e0);

        Edge * e1 = new Edge;
        add_edge(e1);

        Edge * e2 = new Edge;
        add_edge(e2);

        Edge * e3 = new Edge;
        add_edge(e3);

        Edge * e4 = new Edge;
        add_edge(e4);

        Edge * e5 = new Edge;
        add_edge(e5);

        // Save Next/Prev for swapping
        Edge * this_next = Edges[e->Next];
        Edge * this_prev = Edges[e->Prev];
        Edge * twin_next = Edges[Edges[e->Twin]->Next];
        Edge * twin_prev = Edges[Edges[e->Twin]->Prev];

        // Handle constraint curves
        e1->Orientation = e->Orientation;
        e4->Orientation = !e->Orientation;
        Edges[e->Twin]->Orientation = !e->Orientation;
        if (e->Orientation) {
            e1->ConstraintCurve = c;
        } else {
            e1->ConstraintCurve = e->ConstraintCurve;
            e->ConstraintCurve = c;
        }
        Edges[e->Twin]->ConstraintCurve = e1->ConstraintCurve;
        e4->ConstraintCurve = e->ConstraintCurve;

        // Construct Edges
        e0->Node = Points.size() - 1;
        e0->Next = e->Prev;
        e0->Prev = e->Self;
        e0->Twin = e2->Self;
        e0->Mark = false;

        e1->Node = Points.size() - 1;
        e1->Next = e->Next;
        e1->Prev = e2->Self;
        e1->Twin = e->Twin;
        e1->Mark = false;

        e2->Node = Edges[e->Prev]->Node;
        e2->Next = e1->Self;
        e2->Prev = e->Next;
        e2->Twin = e0->Self;
        e2->Mark = false;

        e3->Node = Points.size() - 1;
        e3->Next = Edges[e->Twin]->Prev;
        e3->Prev = e->Twin;
        e3->Twin = e5->Self;
        e3->Mark = false;

        e4->Node = Points.size() - 1;
        e4->Next = Edges[e->Twin]->Next;
        e4->Prev = e5->Self;
        e4->Twin = e->Self;
        e4->Mark = false;

        e5->Node = Edges[Edges[e->Twin]->Prev]->Node;
        e5->Next = e4->Self;
        e5->Prev = Edges[e->Twin]->Next;
        e5->Twin = e3->Self;
        e5->Mark = false;

        Edges[Edges[e->Twin]->Next]->Next = e5->Self;
        Edges[Edges[e->Twin]->Next]->Prev = e4->Self;
        Edges[Edges[e->Twin]->Next]->Mark = false;

        Edges[Edges[e->Twin]->Prev]->Prev = e3->Self;
        Edges[Edges[e->Twin]->Prev]->Mark = false;

        Edges[e->Twin]->Next = e3->Self;
        Edges[e->Twin]->Twin = e1->Self;
        Edges[e->Twin]->Mark = false;

        Edges[e->Next]->Next = e2->Self;
        Edges[e->Next]->Prev = e1->Self;
        Edges[e->Twin]->Mark = false;

        Edges[e->Prev]->Prev = e0->Self;
        Edges[e->Prev]->Mark = false;

        e->Next = e0->Self;
        e->Twin = e4->Self;
        e->Mark = false;

        // Recursive swap
        recursive_swap(e0);
        recursive_swap(e3);
        recursive_swap(this_next);
        recursive_swap(this_prev);
        recursive_swap(twin_next);
        recursive_swap(twin_prev);
    }

    return InsertPointResult::Midpoint;
}

InsertPointResult Mesh::insert_point(Point const vc, Edge *tri) {
    // Find triangle containing point
    LocateTriangleResult result = locate_triangle(vc, tri);

    // Do not insert point if it is duplicate
    if (result == LocateTriangleResult::Point) {
        return InsertPointResult::Duplicate;
    }

    // Test edges in current and adjacent triangles for encroachment
    // These are the only possible edges that are encroached due to empty circumcircle property
    std::vector<Edge *> test_edges;
    test_edges.reserve(9);

    test_edges.push_back(tri);
    test_edges.push_back(Edges[tri->Next]);
    test_edges.push_back(Edges[tri->Prev]);
    for (size_t i = 0; i != 3; ++i) {
        if (test_edges[i]->Self != test_edges[i]->Twin) {
            test_edges.push_back(Edges[Edges[test_edges[i]->Twin]->Next]);
            test_edges.push_back(Edges[Edges[test_edges[i]->Twin]->Prev]);
        }
    }

    // Split all encroached edges
    bool encroached = false;
    for (auto e : test_edges) {
        if (is_encroached(e, vc)) {
            insert_midpoint(e);
            encroached = true;
        }
    }

    // If none are encroached, insert circumcenter
    if (encroached) {
        return InsertPointResult::Midpoint;
    } else if (result == LocateTriangleResult::Interior) {
        // Insert circumcenter
        Edge * next = Edges[tri->Next];
        Edge * prev = Edges[tri->Prev];

        size_t vt = tri->Node;
        size_t vn = next->Node;
        size_t vp = prev->Node;

        Edge *e0 = new Edge;
        add_edge(e0);
        Edge *e1 = new Edge;
        add_edge(e1);
        Edge *e2 = new Edge;
        add_edge(e2);
        Edge *e3 = new Edge;
        add_edge(e3);
        Edge *e4 = new Edge;
        add_edge(e4);
        Edge *e5 = new Edge;
        add_edge(e5);

        Points.push_back(vc);

        e0->Node = Points.size() - 1;
        e0->Next = tri->Self;
        e0->Prev = e1->Self;
        e0->Twin = e5->Self;
        e0->Mark = false;

        e1->Node = vn;
        e1->Next = e0->Self;
        e1->Prev = tri->Self;
        e1->Twin = e2->Self;
        e1->Mark = false;

        e2->Node = Points.size() - 1;
        e2->Next = next->Self;
        e2->Prev = e3->Self;
        e2->Twin = e1->Self;
        e2->Mark = false;

        e3->Node = vp;
        e3->Next = e2->Self;
        e3->Prev = next->Self;
        e3->Twin = e4->Self;
        e3->Mark = false;

        e4->Node = Points.size() - 1;
        e4->Next = prev->Self;
        e4->Prev = e5->Self;
        e4->Twin = e3->Self;
        e4->Mark = false;

        e5->Node = vt;
        e5->Next = e4->Self;
        e5->Prev = prev->Self;
        e5->Twin = e0->Self;
        e5->Mark = false;

        next->Next = e3->Self;
        next->Prev = e2->Self;
        next->Mark = false;

        prev->Next = e5->Self;
        prev->Prev = e4->Self;
        prev->Mark = false;

        tri->Next = e1->Self;
        tri->Prev = e0->Self;
        tri->Mark = false;

        recursive_swap(tri);
        recursive_swap(next);
        recursive_swap(prev);

        return InsertPointResult::Success;
    } else {
        throw std::exception(); // TODO: No triangle could be found containing circumcenter and no boundary edge was encroached by circumcenter
    }
}

Point Mesh::circumcenter(Edge const *e) const {
    Point const node = Points[e->Node];
    Point const prev = Points[Edges[e->Prev]->Node];
    Point const next = Points[Edges[e->Next]->Node];

    double xa = node.X;
    double ya = node.Y;
    double xb = prev.X - xa;
    double yb = prev.Y - ya;
    double xc = next.X - xa;
    double yc = next.Y - ya;

    double d = 2.0 * (xb * yc - yb * xc);
    double db = xb * xb + yb * yb;
    double dc = xc * xc + yc * yc;
    double x = (yc * db - yb * dc) / d + xa;
    double y = (xb * dc - xc * db) / d + ya;

    return Point{x, y};
}