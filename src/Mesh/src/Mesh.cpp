#include "Mesh.hpp"

Mesh::Mesh(Sketch &sketch) {
    Boundary = sketch.boundary();

    for (size_t i = 0; i != sketch.size_curves(); ++i) {
        auto c = sketch.curve(i);
        if (!(c->ForConstruction)) {
            Curves.push_back(c);
        }
    }

    for (size_t i = 0; i != sketch.size_contours(); ++i) {
        Contours.push_back(sketch.contour(i));
    }

    Constraints.push_back(DartConstraint());
}

bool Mesh::are_intersecting(size_t ei, size_t ej) const {
    // TODO, Make more detailed return type enumeration
    if (is_constrained(ei) && constraint_curve(ei) == constraint_curve(ej)) {
        return false;
    }

    Point const v00 = base(ei);
    Point const v01 = tip(ei);
    Point const v10 = base(ej);
    Point const v11 = tip(ej);

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
            dmin = std::fmin(dmin, dx * dx + dy * dy);

            dx = xs0 - xd0 * s - xs1 + xd1 * s;
            dy = ys0 - yd0 * s - ys1 + yd1 * s;
            dmin = std::fmin(dmin, dx * dx + dy * dy);
        }

        s = ((xd0 + xd1) * (xs0 - xs1) + (yd0 + yd1) * (ys0 - ys1)) /
            ((xd0 + xd1) * (xd0 + xd1) + (yd0 + yd1) * (yd0 + yd1));
        if (abs(s) < 1.0 - FLT_EPSILON) {
            dx = xs0 + xd0 * s - xs1 + xd1 * s;
            dy = ys0 + yd0 * s - ys1 + yd1 * s;
            dmin = std::fmin(dmin, dx * dx + dy * dy);

            dx = xs0 - xd0 * s - xs1 - xd1 * s;
            dy = ys0 - yd0 * s - ys1 - yd1 * s;
            dmin = std::fmin(dmin, dx * dx + dy * dy);
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

    for (size_t e = 0; e != Edges.size(); ++e) {
        if (e != prev(next(e))) {
            result = false;
            break;
        }
        if (e != next(prev(e))) {
            result = false;
            break;
        }
        if (e != twin(twin(e))) {
            result = false;
            break;
        }

        if ((e != twin(e))) {
            if (node(e) != node(next(twin(e)))) {
                result = false;
                break;
            }
            if (constraint_curve(e) != constraint_curve(twin(e))) {
                result = false;
                break;
            }
            if (is_constrained(e)) {
                if (orientation(e) == orientation(twin(e))) {
                    result = false;
                    break;
                }
            }

            if (node(e) == node(twin(e))) {
                result = false;
                break;
            }
        }

        if (is_constrained(e)) {
            DartConstraint dc = Constraints[Edges[e].Constraint];
            double tol = length(e) * FLT_EPSILON;
            if (orientation(e)) {
                Point p0 = base(e);
                Point p1 = dc.ConstraintCurve->point(dc.S0);
                if (dist(p0,p1) > tol) {
                    result = false;
                    break;
                }

                p0 = tip(e);
                p1 = dc.ConstraintCurve->point(dc.S1);
                if (dist(p0,p1) > tol) {
                    result = false;
                    break;
                }
            } else {
                Point p0 = base(e);
                Point p1 = dc.ConstraintCurve->point(dc.S1);
                if (dist(p0,p1) > tol) {
                    result = false;
                    break;
                }

                p0 = tip(e);
                p1 = dc.ConstraintCurve->point(dc.S0);
                if (dist(p0,p1) > tol) {
                    result = false;
                    break;
                }
            }
        }
    }

    return result;
}

bool Mesh::find_attached(Point const p, size_t &e_out) {
    double tol = length(e_out) * FLT_EPSILON;

    if (dist(tip(e_out),p) < tol) {
        return true;
    }

    size_t e_in = e_out;

    if (e_out != twin(e_out)) {
        e_out = next(twin(e_out));
        while (e_out != e_in) {
            if (dist(tip(e_out),p) < tol) {
                return true;
            } else if (e_out != twin(e_out)) {
                e_out = next(twin(e_out));
            } else {
                break;
            }
        }
    }

    if (e_out == twin(e_out)) {
        e_out = prev(e_in);
        while (e_out != twin(e_out)) {
            if (dist(base(e_out),p) < tol) {
                e_out = twin(e_out);
                return true;
            } else {
                e_out = prev(twin(e_out));
            }
        }
    }

    return false;
}

bool Mesh::refine() {
    std::vector<double> radii;
    std::vector<double> quality;
    std::vector<size_t> index;

    // TODO: Loop until quality is satisfied
    // TODO: Iteratively decrease the minimum element size until quality is satisfied

    element_quality(radii, quality);
    sort_permutation_ascending(quality, index);
    size_t N = Triangles.size();

    refine_once(index, radii, quality);
    size_t M = Triangles.size();

    while (M > N) {
        N = M;
        element_quality(radii, quality);
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

    element_quality(radii, quality);
    sort_permutation_ascending(quality, index);
    //sort_permutation_descending(radii, index);
    refine_once(index, radii, quality);

    return edges_are_valid();
}

bool Mesh::in_triangle(Point const p, size_t ei) const {
    double xp = p.X;
    double yp = p.Y;

    Point const p0 = base(ei);
    Point const p1 = base(next(ei));
    Point const p2 = base(prev(ei));

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

bool Mesh::is_encroached(Point const p, size_t ei) const {
    /*
        A constrained edge is encroached if a triangle and it's circumcenter lie on opposite sides of the edge.
        This is equivalent to a node being in the diameteral ball of the edge?
        This only occurs if one of the triangles attached to the edge encroaches the edge?
        This only occurs if the angle of the triangles attached to the edge has an angle opposite the edge of greater than 90 degrees.
        Using the dot product, this requires that the normalized dot product < cos(90) = 0
    */

    if (!is_constrained(ei)) {
        return false;
    } else {
        Point const p0 = base(ei);
        Point const p1 = tip(ei);

        double dx0 = p0.X - p.X;
        double dy0 = p0.Y - p.Y;

        double dx1 = p1.X - p.X;
        double dy1 = p1.Y - p.Y;

        double dot = dx0 * dx1 + dy0 * dy1;
        double tol = std::sqrt(dx0 * dx0 + dy0 * dy0) * std::sqrt(dx1 * dx1 + dy1 * dy1) * FLT_EPSILON;

        return (dot < tol);
    }
}

bool Mesh::is_locally_optimal(size_t ei) const {
    /*
        See Chapter 3.7 of "Triangulations and Applications" by Øyvind Hjelle and Morten Dæhlen
    */

    if (is_constrained(ei)) {
        return true;
    } else {
        Point const p3 = base(ei);
        Point const p2 = base(prev(ei));
        Point const p1 = base(twin(ei));
        Point const p4 = base(prev(twin(ei)));

        double v1x = p3.X - p2.X;
        double v1y = p3.Y - p2.Y;
        double v2x = p1.X - p2.X;
        double v2y = p1.Y - p2.Y;
        double v3x = p1.X - p4.X;
        double v3y = p1.Y - p4.Y;
        double v4x = p3.X - p4.X;
        double v4y = p3.Y - p4.Y;

        double d1 = std::sqrt(v1x * v1x + v1y * v1y);
        double d2 = std::sqrt(v2x * v2x + v2y * v2y);
        double d3 = std::sqrt(v3x * v3x + v3y * v3y);
        double d4 = std::sqrt(v4x * v4x + v4y * v4y);
        double tol = -(d1 * d2 * d3 * d4 * FLT_EPSILON);

        double sina = v1x * v2y - v1y * v2x;
        double sinb = v3x * v4y - v3y * v4x;
        double cosa = v1x * v2x + v1y * v2y;
        double cosb = v3x * v4x + v3y * v4y;
        double cct = sina * cosb + cosa * sinb;

        return cct >= 0; // cct >= 0 > tol?
    }
}

bool Mesh::is_protruding(size_t ei) const {
    //
    //    See Chapter 1.4 of "Triangulations and Applications" by Øyvind Hjelle and Morten Dæhlen
    //

    Point const p0 = base(prev(ei));
    Point const p1 = base(ei);
    Point const p2 = base(next(ei));

    double v1x = p2.X - p1.X;
    double v1y = p2.Y - p1.Y;

    double v0x = p0.X - p1.X;
    double v0y = p0.Y - p1.Y;

    double area = v1x * v0y - v1y * v0x;

    if (area > 0.0) {
        // Make sure no boundary points interior to triangle
        // Calculate barrycentric coordinates of p

        size_t nxt = next(next(ei));
        while (nxt != prev(ei)) {
            Point const p4 = point(nxt);

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
                nxt = next(nxt);
            }
        }
        return true;
    } else {
        return false;
    }
}

bool Mesh::is_valid(size_t ei) const {
    bool value = true;

    value = value && (ei == prev(next(ei)));
    value = value && (ei == next(prev(ei)));
    value = value && (ei == twin(twin(ei)));

    return value;
}

bool Mesh::recursive_swap(size_t ei) {
    // TODO, May need to have two different recursive swap methods, one for midpoint insertion and one for circumcenter insertion
    if (!is_locally_optimal(ei) && swap(ei)) {
        size_t enext = next(ei);
        size_t eprev = prev(ei);
        size_t tnext = next(twin(ei));
        size_t tprev = prev(twin(ei));

        recursive_swap(enext);
        recursive_swap(eprev);
        recursive_swap(tnext);
        recursive_swap(tprev);

        return true;
    } else {
        return false;
    }
}

bool Mesh::swap(size_t ei) {
    if (!is_constrained(ei)) {
        Edge &e0 = Edges[ei];
        Edge &e1 = Edges[e0.Next];
        Edge &e2 = Edges[e0.Prev];
        Edge &e5 = Edges[e0.Twin];
        Edge &e3 = Edges[e5.Next];
        Edge &e4 = Edges[e5.Prev];

        e0.Node = e2.Node;
        e0.Next = e4.Self;
        e0.Prev = e1.Self;
        e0.Mark = false;

        e5.Node = e4.Node;
        e5.Next = e2.Self;
        e5.Prev = e3.Self;
        e5.Mark = false;

        e1.Next = e0.Self;
        e1.Prev = e4.Self;
        e1.Mark = false;

        e2.Next = e3.Self;
        e2.Prev = e0.Twin;
        e2.Mark = false;

        e3.Next = e0.Twin;
        e3.Prev = e2.Self;
        e3.Mark = false;

        e4.Next = e1.Self;
        e4.Prev = e0.Self;
        e4.Mark = false;

        return true;
    } else {
        Edges[ei].Mark = false;
        return false;
    }
}

double Mesh::circumradius(size_t ei) const {
    Point const slf = base(ei);
    Point const prv = base(prev(ei));
    Point const nxt = base(next(ei));

    double xa = slf.X;
    double ya = slf.Y;
    double xb = prv.X - xa;
    double yb = prv.Y - ya;
    double xc = nxt.X - xa;
    double yc = nxt.Y - ya;
    xa = xb - xc;
    ya = yb - yc;

    double den = 2.0 * abs(xb * yc - yb * xc);
    double num = std::sqrt(xa * xa + ya * ya) * std::sqrt(xb * xb + yb * yb) * std::sqrt(xc * xc + yc * yc);

    return num / den;
}

double Mesh::length(size_t ei) const {
    Point const p0 = base(ei);
    Point const p1 = tip(ei);

    double dx = p0.X - p1.X;
    double dy = p0.Y - p1.Y;

    return std::sqrt(dx * dx + dy * dy);
}

double Mesh::shortest_edge_length(size_t ei) const { // TODO: Rename Edge class to Dart
    Point const slf = base(ei);
    Point const prv = base(prev(ei));
    Point const nxt = base(next(ei));

    double x0 = slf.X;
    double y0 = slf.Y;
    double x1 = nxt.X;
    double y1 = nxt.Y;
    double x2 = prv.X;
    double y2 = prv.Y;

    double dx = x0 - x1;
    double dy = y0 - y1;
    double dl = dx * dx + dy * dy;

    dx = x1 - x2;
    dy = y1 - y2;
    dl = std::fmin(dl, dx * dx + dy * dy);

    dx = x2 - x0;
    dy = y2 - y0;
    dl = std::fmin(dl, dx * dx + dy * dy);

    return std::sqrt(dl);
}

size_t Mesh::num_edges() const {
    size_t count = 0;
    for (size_t e = 0; e != Edges.size(); ++e) {
        count += (e == twin(e) ? 2 : 1);
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
        Curve const *cc = Boundary->curve(i); //clone() to prevent alteration of input Contour when Edge is split
        bool dir = Boundary->orientation(i);
        Edge &e = new_edge(Points.size(), Constraints.size(), dir); // TODO: change to new_edge(Edge)
        e.Twin = e.Self;

        Constraints.push_back(DartConstraint(0.0, 1.0, cc));
        Points.push_back(Point(dir ? cc->start() : cc->end()));
    }

    // Set Next/Prev edges from Boundary
    size_t Ne = 0;
    for (size_t i = 0; i != Edges.size(); ++i) {
        size_t j = (i + 1) % Edges.size();
        Edges[i].Next = Edges[j].Self;
        Edges[j].Prev = Edges[i].Self;
    }

    // TODO: Length Based Refinement of Boundary curves

    // Some edges may intersect due to discretization error
    // If two edges intersect, split the longest edge
    bool any_split = true;
    while (any_split) {
        any_split = false;
        for (size_t i = 0; i != Edges.size(); ++i) {
            for (size_t j = i + 1; j != Edges.size(); ++j) {
                if (are_intersecting(i, j)) {
                    any_split = true;
                    if (length(i) > length(j)) {
                        split_edge(i);
                    } else {
                        split_edge(j);
                    }
                }
            }
        }
    }
}

void Mesh::element_quality(std::vector<double> &radii, std::vector<double> &quality) {
    radii.resize(0);
    quality.resize(0);

    radii.reserve(Triangles.size());
    quality.reserve(Triangles.size());
    for (size_t i = 0; i < Triangles.size(); ++i) {
        double r = circumradius(Triangles[i]);
        double l = shortest_edge_length(Triangles[i]);

        radii.push_back(r);
        quality.push_back(l / r);
    }
}

void Mesh::get_triangles() {
    Triangles.resize(0);
    Triangles.reserve(2 * num_points());

    mark_triangles();
    for (auto &e : Edges) {
        if (e.Mark) {
            Triangles.push_back(e.Self);
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
    std::vector<Curve const *> interior;
    for (auto c : Curves) {
        bool on_exterior = false;
        for (size_t i = 0; i != Boundary->size(); ++i) {
            if (c == Boundary->curve(i)) {
                on_exterior = true;
                break;
            }
        }

        if (!on_exterior) {
            interior.push_back(c);
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
    std::vector<DartConstraint> queue;
    for (size_t i = 0; i != interior.size(); ++i) {
        queue.push_back(DartConstraint(0.0, 1.0, interior[i]));
        while (queue.size() != 0) {
            DartConstraint &dc = queue.back();
            Point p0 = dc.ConstraintCurve->point(dc.S0); // TODO: write Point Curve::point(double) and differentiate from Vertex Curve::vertex(double)
            Point p1 = dc.ConstraintCurve->point(dc.S1);

            size_t ei = Edges.size() - 1;
            LocateTriangleResult result = locate_triangle(p0, ei);

            if (result != LocateTriangleResult::Point) {
                throw std::exception();
            }

            if (find_attached(p1, ei)) {
                Edge &e = Edges[ei];
                e.Constraint = Constraints.size();
                e.Orientation = true;

                Edge &et = Edges[e.Twin];
                et.Constraint = Constraints.size();
                et.Orientation = false;

                Constraints.push_back(queue.back());
                queue.pop_back();
            } else {
                double s0 = dc.S0;
                double s1 = dc.S1;
                Curve const *cc = dc.ConstraintCurve;

                double sn = (s0 + s1) / 2.0;
                Point const p = cc->point(sn);

                dc.S1 = sn;
                queue.push_back(DartConstraint(sn, s1, cc));

                while (insert_point(p) == InsertPointResult::Midpoint);
            }
        }
    }

    split_encroached_edges();
}

void Mesh::mark_triangles() {
    for (size_t i = 0; i != Edges.size(); ++i) {
        Edges[i].Mark = true;
    }

    for (size_t i = 0; i != Edges.size(); ++i) {
        if (Edges[i].Mark && Edges[next(i)].Mark && Edges[prev(i)].Mark) {
            Edges[next(i)].Mark = false;
            Edges[prev(i)].Mark = false;
        }
    }
}

void Mesh::refine_once(std::vector<size_t> index, std::vector<double> radii, std::vector<double> quality) {
    for (size_t i = 0; i != Triangles.size(); ++i) {
        size_t j = index[i];
        if ((triangle(j).Mark) && ((radii[j] > MaximumElementSize) || (radii[j] > MinimumElementSize && quality[j] < MinimumElementQuality))) {
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

    for (size_t e = 0; e != Edges.size(); ++e) {
        Point const v0 = base(e);
        Point const v1 = base(next(e));
        Point const v2 = base(next(next(e)));
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

void Mesh::split_edge(size_t ei) {
    /*
        Splits edge into two edges at the midpoint without creating any new triangles.
        Used for initial polygon refinement.
    */

    // TODO: Refactor into split_constrained_edge
    size_t c{0};
    if (is_constrained(ei)) {
        c = Constraints.size();

        DartConstraint &dc = Constraints[Edges[ei].Constraint];
        double s0 = dc.S0;
        double s1 = dc.S1;
        Curve const *cc = dc.ConstraintCurve;

        double sn = (s0 + s1) / 2.0;
        dc.S1 = sn;

        Constraints.push_back(DartConstraint(sn, s1, cc));
        Points.push_back(cc->point(sn));
    } else { // Unconstrained Edge
        Point const p0 = base(ei);
        Point const p1 = tip(ei);

        Points.push_back(Point((p0.X + p1.X) / 2.0, (p0.Y + p1.Y) / 2.0));
    }

    if (ei == twin(ei)) { // Boundary Edge, TODO: write bool Edge::is_boundary()
        size_t itr = new_edges(1);
        Edge &newe = Edges[--itr];
        Edge &e = Edges[ei];
        Edge &nxt = Edges[e.Next];

        // Constraint Curve
        newe.Orientation = e.Orientation;
        if (e.Orientation) {
            newe.Constraint = c;
        } else {
            newe.Constraint = e.Constraint;
            e.Constraint = c;
        }

        // Connectivity
        newe.Node = Points.size() - 1;
        newe.Next = e.Next;
        newe.Prev = e.Self;
        newe.Twin = newe.Self;
        newe.Mark = false;

        nxt.Prev = newe.Self;
        nxt.Mark = false;

        e.Next = newe.Self;
        e.Mark = false;
    } else {
        throw std::exception(); // Function should only be called on constrained edges
    }
}

void Mesh::split_encroached_edges() {
    bool any_split = true;
    while (any_split) {
        any_split = false;
        for (size_t i = 0; i != Edges.size(); ++i) {
            if (is_constrained(i)) {
                if (is_encroached(base(prev(i)), i)) {
                    any_split = true;
                    insert_midpoint(i);
                }
            }
        }
    }
}

void Mesh::triangulate_boundary_polygon() {
    Edges.reserve(3 * num_points());

    size_t i = 0;
    while (i != next(next(next(i)))) {
        if (is_protruding(i)) {
            size_t itr = new_edges(2);
            Edge &e1 = Edges[--itr];
            Edge &e0 = Edges[--itr];

            Edge &ei = Edges[i];
            Edge &nxt = Edges[ei.Next];
            Edge &prv = Edges[ei.Prev];

            Edge &prvprv = Edges[prv.Prev];

            // Edge of new triangle
            e0.Node = nxt.Node;
            e0.Prev = ei.Self;
            e0.Next = ei.Prev;
            e0.Twin = e1.Self;

            // Twin edge, part of new polygonal boundary
            e1.Node = prv.Node;
            e1.Next = ei.Next;
            e1.Prev = prv.Prev;
            e1.Twin = e0.Self;

            // Update polygonal boundary
            nxt.Prev = e1.Self;
            ei.Next = e0.Self;
            prvprv.Next = e1.Self;
            prv.Prev = e0.Self;

            // Next edge
            i = next(e1.Self);
        } else {
            i = next(i);
        }
    }

    // Edge swap to make triangulation Delaunay
    for (size_t i = 0; i != Edges.size(); ++i) {
        recursive_swap(i);
    }

    split_encroached_edges();
}

LocateTriangleResult Mesh::locate_triangle(Point const p, size_t &ei) const {
    double xp = p.X;
    double yp = p.Y;

    Point p0 = base(ei);
    Point p1 = base(next(ei));
    Point p2 = base(prev(ei));

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

    double tol_a = FLT_EPSILON * (dx20 * dy01 - dy20 * dx01);
    double tol_l = FLT_EPSILON * sqrt(dx20 * dy01 - dy20 * dx01);

    double dist0 = sqrt(dx0p * dx0p + dy0p * dy0p);
    double dist1 = sqrt(dx1p * dx1p + dy1p * dy1p);
    double dist2 = sqrt(dx2p * dx2p + dy2p * dy2p);

    double area01 = dx0p * dy1p - dx1p * dy0p;
    double area12 = dx1p * dy2p - dx2p * dy1p;
    double area20 = dx2p * dy0p - dx0p * dy2p;

    if (dist0 < tol_l) {
        return LocateTriangleResult::Point;
    } else if (dist1 < tol_l) {
        ei = next(ei);
        return LocateTriangleResult::Point;
    } else if (dist2 < tol_l) {
        ei = prev(ei);
        return LocateTriangleResult::Point;
    } else  if (area01 > tol_a && area12 > tol_a && area20 > tol_a) {
        return LocateTriangleResult::Interior;
    } else if (area01 < -tol_a && ei != twin(ei)) {
        ei = twin(ei);
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
    } else if (area12 < -tol_a && next(ei) != twin(next(ei))) {
        ei = twin(next(ei));
        p0 = p2;

        dx0p = dx2p;
        dy0p = dy2p;

        dx01 = -dx12;
        dy01 = -dy12;
    } else if (area20 < -tol_a && prev(ei) != twin(prev(ei))) {
        ei = twin(prev(ei));
        p1 = p2;

        dx1p = dx2p;
        dy1p = dy2p;

        dx01 = -dx20;
        dy01 = -dy20;
    } else if (area01 > -tol_a && area12 > tol_a && area20 > tol_a) {
        ei = twin(ei);
        return LocateTriangleResult::Interior;
    } else if (area01 > tol_a && area12 > -tol_a && area20 > tol_a) {
        ei = twin(next(ei));
        return LocateTriangleResult::Interior;
    } else if (area01 > tol_a && area12 > tol_a && area20 > -tol_a) {
        ei = twin(prev(ei));
        return LocateTriangleResult::Interior;
    } else if (area01 < -tol_a) {
        return LocateTriangleResult::Exterior;
    } else if (area12 < -tol_a) {
        ei = next(ei);
        return LocateTriangleResult::Exterior;
    } else if (area20 < -tol_a) {
        ei = prev(ei);
        return LocateTriangleResult::Exterior;
    } else {
        throw std::exception();
    }

    while (true) { // After first iteration, area01 > 0
        p2 = base(prev(ei));

        dx2p = p2.X - xp;
        dy2p = p2.Y - yp;

        dx12 = p1.X - p2.X;
        dy12 = p1.Y - p2.Y;

        dx20 = p2.X - p0.X;
        dy20 = p2.Y - p0.Y;

        tol_a = FLT_EPSILON * (dx20 * dy01 - dy20 * dx01);
        tol_l = FLT_EPSILON * sqrt(dx20 * dy01 - dy20 * dx01);

        dist2 = sqrt(dx2p * dx2p + dy2p * dy2p);

        area12 = dx1p * dy2p - dx2p * dy1p;
        area20 = dx2p * dy0p - dx0p * dy2p;

        if (dist2 < tol_l) {
            ei = prev(ei);
            return LocateTriangleResult::Point;
        } else if (area12 > tol_a && area20 > tol_a) {
            return LocateTriangleResult::Interior;
        } else if (area12 < -tol_a && next(ei) != twin(next(ei))) {
            ei = twin(next(ei));
            p0 = p2;

            dx0p = dx2p;
            dy0p = dy2p;

            dx01 = -dx12;
            dy01 = -dy12;
            continue;
        } else if (area20 < -tol_a && prev(ei) != twin(prev(ei))) {
            ei = twin(prev(ei));
            p1 = p2;

            dx1p = dx2p;
            dy1p = dy2p;

            dx01 = -dx20;
            dy01 = -dy20;
            continue;
        } else if (area12 > -tol_a && area20 > tol_a) {
            ei = twin(next(ei));
            return LocateTriangleResult::Interior;
        } else if (area12 > tol_a && area20 > -tol_a) {
            ei = twin(prev(ei));
            return LocateTriangleResult::Interior;
        } else if (area12 < -tol_a) {
            ei = next(ei);
            return LocateTriangleResult::Exterior;
        } else if (area20 < -tol_a) {
            ei = prev(ei);
            return LocateTriangleResult::Exterior;
        } else {
            throw std::exception();
        }
    }
}

InsertPointResult Mesh::insert_circumcenter(size_t ei) {
    return insert_point(circumcenter(ei), ei);
}

InsertPointResult Mesh::insert_midpoint(size_t ei) {
    /*
        Splits edge into two edges and creates two new triangles.
    */

    size_t c{0};
    if (is_constrained(ei)) { // Constrained Edge
        c = Constraints.size();

        DartConstraint &dc = Constraints[Edges[ei].Constraint];
        double s0 = dc.S0;
        double s1 = dc.S1;
        Curve const *cc = dc.ConstraintCurve;

        double sn = (s0 + s1) / 2.0;
        dc.S1 = sn;

        Constraints.push_back(DartConstraint(sn, s1, cc));
        Points.push_back(cc->point(sn));
    } else { // Unconstrained Edge
        Point const p0 = base(ei);
        Point const p1 = tip(ei);
        Points.push_back(Point((p0.X + p1.X) / 2.0, (p0.Y + p1.Y) / 2.0));
    }

    if (ei == twin(ei)) { // Boundary Edge
        size_t itr = new_edges(3);
        Edge &e2 = Edges[--itr];
        Edge &e1 = Edges[--itr];
        Edge &e0 = Edges[--itr];

        Edge &e = Edges[ei];
        Edge &nxt = Edges[e.Next];
        Edge &prv = Edges[e.Prev];

        // Handle constraint curves
        e2.Orientation = e.Orientation;
        if (e.Orientation) {
            e2.Constraint = c;
        } else {
            e2.Constraint = e.Constraint;
            e.Constraint = c;
        }

        // Construct edges
        e0.Node = prv.Node;
        e0.Next = e2.Self;
        e0.Prev = e.Next;
        e0.Twin = e1.Self;
        e0.Mark = false;

        e1.Node = Points.size() - 1; // TODO: write new_point(Point) method
        e1.Next = e.Prev;
        e1.Prev = e.Self;
        e1.Twin = e0.Self;
        e1.Mark = false;

        e2.Node = Points.size() - 1;
        e2.Next = e.Next;
        e2.Prev = e0.Self;
        e2.Twin = e2.Self;
        e2.Mark = false;

        nxt.Next = e0.Self;
        nxt.Prev = e2.Self;
        nxt.Mark = false;

        prv.Prev = e1.Self;
        prv.Mark = false;

        e.Next = e1.Self;
        e.Mark = false;

        // Recursive swap
        recursive_swap(e0.Self);
        recursive_swap(nxt.Self);
        recursive_swap(prv.Self);
    } else { // Interior Edge
        size_t itr = new_edges(6);
        Edge &e5 = Edges[--itr];
        Edge &e4 = Edges[--itr];
        Edge &e3 = Edges[--itr];
        Edge &e2 = Edges[--itr];
        Edge &e1 = Edges[--itr];
        Edge &e0 = Edges[--itr];

        Edge &e = Edges[ei];
        Edge &nxt = Edges[e.Next];
        Edge &prv = Edges[e.Prev];
        Edge &twn = Edges[e.Twin];
        Edge &tnxt = Edges[twn.Next];
        Edge &tprv = Edges[twn.Prev];

        // Handle constraint curves
        e1.Orientation = e.Orientation;
        e4.Orientation = !e.Orientation;
        twn.Orientation = !e.Orientation;
        if (e.Orientation) {
            e1.Constraint = c;
        } else {
            e1.Constraint = e.Constraint;
            e.Constraint = c;
        }
        twn.Constraint = e1.Constraint;
        e4.Constraint = e.Constraint;

        // Construct Edges
        e0.Node = Points.size() - 1;
        e0.Next = e.Prev;
        e0.Prev = e.Self;
        e0.Twin = e2.Self;
        e0.Mark = false;

        e1.Node = Points.size() - 1;
        e1.Next = e.Next;
        e1.Prev = e2.Self;
        e1.Twin = e.Twin;
        e1.Mark = false;

        e2.Node = prv.Node;
        e2.Next = e1.Self;
        e2.Prev = e.Next;
        e2.Twin = e0.Self;
        e2.Mark = false;

        e3.Node = Points.size() - 1;
        e3.Next = twn.Prev;
        e3.Prev = e.Twin;
        e3.Twin = e5.Self;
        e3.Mark = false;

        e4.Node = Points.size() - 1;
        e4.Next = twn.Next;
        e4.Prev = e5.Self;
        e4.Twin = e.Self;
        e4.Mark = false;

        e5.Node = tprv.Node;
        e5.Next = e4.Self;
        e5.Prev = twn.Next;
        e5.Twin = e3.Self;
        e5.Mark = false;

        tnxt.Next = e5.Self;
        tnxt.Prev = e4.Self;
        tnxt.Mark = false;

        tprv.Prev = e3.Self;
        tprv.Mark = false;

        twn.Next = e3.Self;
        twn.Twin = e1.Self;
        twn.Mark = false;

        nxt.Next = e2.Self;
        nxt.Prev = e1.Self;
        nxt.Mark = false;

        prv.Prev = e0.Self;
        prv.Mark = false;

        e.Next = e0.Self;
        e.Twin = e4.Self;
        e.Mark = false;

        // Recursive swap
        recursive_swap(e0.Self);
        recursive_swap(e3.Self);
        recursive_swap(nxt.Self);
        recursive_swap(prv.Self);
        recursive_swap(tnxt.Self);
        recursive_swap(tprv.Self);
    }

    return InsertPointResult::Midpoint;
}

InsertPointResult Mesh::insert_point(Point const vc, size_t ei) {
    // Find triangle containing point
    LocateTriangleResult result = locate_triangle(vc, ei);

    // Do not insert point if it is duplicate
    if (result == LocateTriangleResult::Point) {
        return InsertPointResult::Duplicate;
    }

    // Test edges in current and adjacent triangles for encroachment
    // These are the only possible edges that are encroached due to empty circumcircle property
    std::vector<size_t> test_edges;
    test_edges.reserve(9);

    test_edges.push_back(ei);
    test_edges.push_back(next(ei));
    test_edges.push_back(prev(ei));
    for (size_t i = 0; i != 3; ++i) {
        Edge e = Edges[test_edges[i]];
        if (e.Self != e.Twin) {
            test_edges.push_back(twin(e).Next);
            test_edges.push_back(twin(e).Prev);
        }
    }

    // Split all encroached edges
    bool encroached = false;
    for (auto e : test_edges) {
        if (is_encroached(vc, e)) {
            insert_midpoint(e);
            encroached = true;
        }
    }

    // If none are encroached, insert circumcenter
    if (encroached) {
        return InsertPointResult::Midpoint;
    } else if (result == LocateTriangleResult::Interior) {
        size_t itr = new_edges(6);
        Edge &e5 = Edges[--itr];
        Edge &e4 = Edges[--itr];
        Edge &e3 = Edges[--itr];
        Edge &e2 = Edges[--itr];
        Edge &e1 = Edges[--itr];
        Edge &e0 = Edges[--itr];

        Edge &tri = Edges[ei];
        Edge &nxt = Edges[tri.Next];//next(tri);
        Edge &prv = Edges[tri.Prev];//prev(tri);

        size_t vt = node(tri);
        size_t vn = node(nxt);
        size_t vp = node(prv);

        Points.push_back(vc);

        e0.Node = Points.size() - 1;
        e0.Next = tri.Self;
        e0.Prev = e1.Self;
        e0.Twin = e5.Self;
        e0.Mark = false;

        e1.Node = vn;
        e1.Next = e0.Self;
        e1.Prev = tri.Self;
        e1.Twin = e2.Self;
        e1.Mark = false;

        e2.Node = Points.size() - 1;
        e2.Next = nxt.Self;
        e2.Prev = e3.Self;
        e2.Twin = e1.Self;
        e2.Mark = false;

        e3.Node = vp;
        e3.Next = e2.Self;
        e3.Prev = nxt.Self;
        e3.Twin = e4.Self;
        e3.Mark = false;

        e4.Node = Points.size() - 1;
        e4.Next = prv.Self;
        e4.Prev = e5.Self;
        e4.Twin = e3.Self;
        e4.Mark = false;

        e5.Node = vt;
        e5.Next = e4.Self;
        e5.Prev = prv.Self;
        e5.Twin = e0.Self;
        e5.Mark = false;

        nxt.Next = e3.Self;
        nxt.Prev = e2.Self;
        nxt.Mark = false;

        prv.Next = e5.Self;
        prv.Prev = e4.Self;
        prv.Mark = false;

        tri.Next = e1.Self;
        tri.Prev = e0.Self;
        tri.Mark = false;

        recursive_swap(tri.Self);
        recursive_swap(nxt.Self);
        recursive_swap(prv.Self);

        return InsertPointResult::Success;
    } else {
        throw std::exception(); // TODO: No triangle could be found containing circumcenter and no boundary edge was encroached by circumcenter
    }
}

Point Mesh::circumcenter(size_t ei) const {
    Point const slf = base(ei);
    Point const prv = base(prev(ei));
    Point const nxt = base(next(ei));

    double xa = slf.X;
    double ya = slf.Y;
    double xb = prv.X - xa;
    double yb = prv.Y - ya;
    double xc = nxt.X - xa;
    double yc = nxt.Y - ya;

    double d = 2.0 * (xb * yc - yb * xc);
    double db = xb * xb + yb * yb;
    double dc = xc * xc + yc * yc;
    double x = (yc * db - yb * dc) / d + xa;
    double y = (xb * dc - xc * db) / d + ya;

    return Point{x, y};
}