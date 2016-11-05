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

size_t Mesh::num_edges() const {
    size_t count = 0;
    for (size_t i = 0; i < Edges.size(); ++i) {
        count += (Edges[i] == Edges[i]->Twin ? 2 : 1);
    }
    count /= 2;

    return count;
}

void Mesh::delete_me() {
    for (auto &ed : Edges) {
        delete ed;
    }
    Edges.clear();

    for (auto &p : Points) {
        delete p;
    }
    Points.clear();
}

// Topological Queries
LocateTriangleResult Mesh::locate_triangle(Point const *p, Edge *&e) const {
    Edge const *ec = e;

    LocateTriangleResult result = locate_triangle(p, ec);

    e = const_cast<Edge *>(ec);

    return result;
}

LocateTriangleResult Mesh::locate_triangle(Point const *p, Edge const *&e) const {
    double xp = p->X;
    double yp = p->Y;

    const Point *p0 = e->node();
    const Point *p1 = e->next()->node();
    const Point *p2 = e->prev()->node();

    if (*p == *p0) {
        return LocateTriangleResult::Point;
    } else if (*p == *p1) {
        e = e->next();
        return LocateTriangleResult::Point;
    } else if (*p == *p2) {
        e = e->prev();
        return LocateTriangleResult::Point;
    }

    double dx0p = p0->X - xp;
    double dy0p = p0->Y - yp;

    double dx1p = p1->X - xp;
    double dy1p = p1->Y - yp;

    double dx2p = p2->X - xp;
    double dy2p = p2->Y - yp;

    double dx01 = p0->X - p1->X;
    double dy01 = p0->Y - p1->Y;

    double dx12 = p1->X - p2->X;
    double dy12 = p1->Y - p2->Y;

    double dx20 = p2->X - p0->X;
    double dy20 = p2->Y - p0->Y;

    double tol = FLT_EPSILON * (dx20 * dy01 - dy20 * dx01);

    double area01 = dx0p * dy1p - dx1p * dy0p;
    double area12 = dx1p * dy2p - dx2p * dy1p;
    double area20 = dx2p * dy0p - dx0p * dy2p;

    if (area01 > tol && area12 > tol && area20 > tol) {
        return LocateTriangleResult::Interior;
    } else if (area01 < -tol && e != e->twin()) {
        e = e->twin();
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
    } else if (area12 < -tol && e->next() != e->next()->twin()) {
        e = e->next()->twin();
        p0 = p2;

        dx0p = dx2p;
        dy0p = dy2p;

        dx01 = -dx12;
        dy01 = -dy12;
    } else if (area20 < -tol && e->prev() != e->prev()->twin()) {
        e = e->prev()->twin();
        p1 = p2;

        dx1p = dx2p;
        dy1p = dy2p;

        dx01 = -dx20;
        dy01 = -dy20;
    } else if (area01 > -tol && area12 > tol && area20 > tol) {
        e = e->twin();
        return LocateTriangleResult::Interior;
    } else if (area01 > tol && area12 > -tol && area20 > tol) {
        e = e->next()->twin();
        return LocateTriangleResult::Interior;
    } else if (area01 > tol && area12 > tol && area20 > -tol) {
        e = e->prev()->twin();
        return LocateTriangleResult::Interior;
    } else if (area01 < -tol) {
        return LocateTriangleResult::Exterior;
    } else if (area12 < -tol) {
        e = e->next();
        return LocateTriangleResult::Exterior;
    } else if (area20 < -tol) {
        e = e->prev();
        return LocateTriangleResult::Exterior;
    } else {
        throw std::exception();
    }

    while (true) { // After first iteration, area01 > 0
        p2 = e->prev()->node();

        if (*p == *p2) {
            e = e->prev();
            return LocateTriangleResult::Point;
        }

        dx2p = p2->X - xp;
        dy2p = p2->Y - yp;

        dx12 = p1->X - p2->X;
        dy12 = p1->Y - p2->Y;

        dx20 = p2->X - p0->X;
        dy20 = p2->Y - p0->Y;

        tol = FLT_EPSILON * (dx20 * dy01 - dy20 * dx01);

        area12 = dx1p * dy2p - dx2p * dy1p;
        area20 = dx2p * dy0p - dx0p * dy2p;

        if (area12 > tol && area20 > tol) {
            return LocateTriangleResult::Interior;
        } else if (area12 < -tol && e->next() != e->next()->twin()) {
            e = e->next()->twin();
            p0 = p2;

            dx0p = dx2p;
            dy0p = dy2p;

            dx01 = -dx12;
            dy01 = -dy12;
            continue;
        } else if (area20 < -tol && e->prev() != e->prev()->twin()) {
            e = e->prev()->twin();
            p1 = p2;

            dx1p = dx2p;
            dy1p = dy2p;

            dx01 = -dx20;
            dy01 = -dy20;
            continue;
        } else if (area12 > -tol && area20 > tol) {
            e = e->next()->twin();
            return LocateTriangleResult::Interior;
        } else if (area12 > tol && area20 > -tol) {
            e = e->prev()->twin();
            return LocateTriangleResult::Interior;
        } else if (area12 < -tol) {
            e = e->next();
            return LocateTriangleResult::Exterior;
        } else if (area20 < -tol) {
            e = e->prev();
            return LocateTriangleResult::Exterior;
        } else {
            throw std::exception();
        }
    }
}

bool Mesh::in_triangle(Point const *p, Edge const *&e) const {
    double xp = p->X;
    double yp = p->Y;

    Point const *p0 = e->node();
    Point const *p1 = e->next()->node();
    Point const *p2 = e->prev()->node();

    double dx0p = p0->X - xp;
    double dy0p = p0->Y - yp;

    double dx1p = p1->X - xp;
    double dy1p = p1->Y - yp;

    double dx2p = p2->X - xp;
    double dy2p = p2->Y - yp;

    double dx01 = p0->X - p1->X;
    double dy01 = p0->Y - p1->Y;

    double dx12 = p1->X - p2->X;
    double dy12 = p1->Y - p2->Y;

    double dx20 = p2->X - p0->X;
    double dy20 = p2->Y - p0->Y;

    double area012 = dx01 * dy12 - dy01 * dx12;

    double tol = FLT_EPSILON * area012;

    double area01p = dx0p * dy1p - dx1p * dy0p;
    double area12p = dx1p * dy2p - dx2p * dy1p;
    double area20p = dx2p * dy0p - dx0p * dy2p;

    return (area01p > -tol && area12p > -tol && area20p > -tol);
}

// Point Insertion
InsertPointResult Mesh::insert_point(Point const *vc, Edge *tri) {
    // Find triangle containing point
    LocateTriangleResult result = locate_triangle(vc, tri);

    // Do not insert point if it is duplicate
    if (result == LocateTriangleResult::Point) {
        return InsertPointResult::Duplicate;
    }

    // Test edges in current and adjacent triangles for encroachment
    // These are the only possible edges that are encroached due to empty circumcircle property
    // TODO: write is_encroached method
    std::vector<Edge *> test_edges;
    test_edges.reserve(9);

    test_edges.push_back(tri);
    test_edges.push_back(tri->Next);
    test_edges.push_back(tri->Prev);
    for (size_t i = 0; i != 3; ++i) {
        if (test_edges[i] != test_edges[i]->Twin) {
            test_edges.push_back(test_edges[i]->Twin->Next);
            test_edges.push_back(test_edges[i]->Twin->Prev);
        }
    }

    // Split all encroached edges
    bool encroached = false;
    for (auto e : test_edges) {
        if (e->is_encroached(vc)) {
            insert_midpoint(e);
            encroached = true;
        }
    }

    // If none are encroached, insert circumcenter
    if (encroached) {
        return InsertPointResult::Midpoint;
    } else if (result == LocateTriangleResult::Interior) {
        // Insert circumceneter
        Edge * next = tri->Next;
        Edge * prev = tri->Prev;

        Point const *vt = tri->Node;
        Point const *vn = next->Node;
        Point const *vp = prev->Node;

        Edge *e0 = new Edge;
        Edge *e1 = new Edge;
        Edge *e2 = new Edge;
        Edge *e3 = new Edge;
        Edge *e4 = new Edge;
        Edge *e5 = new Edge;

        Points.push_back(vc);
        Edges.push_back(e0);
        Edges.push_back(e1);
        Edges.push_back(e2);
        Edges.push_back(e3);
        Edges.push_back(e4);
        Edges.push_back(e5);

        e0->Node = vc;
        e0->Next = tri;
        e0->Prev = e1;
        e0->Twin = e5;
        e0->Mark = false;

        e1->Node = vn;
        e1->Next = e0;
        e1->Prev = tri;
        e1->Twin = e2;
        e1->Mark = false;

        e2->Node = vc;
        e2->Next = next;
        e2->Prev = e3;
        e2->Twin = e1;
        e2->Mark = false;

        e3->Node = vp;
        e3->Next = e2;
        e3->Prev = next;
        e3->Twin = e4;
        e3->Mark = false;

        e4->Node = vc;
        e4->Next = prev;
        e4->Prev = e5;
        e4->Twin = e3;
        e4->Mark = false;

        e5->Node = vt;
        e5->Next = e4;
        e5->Prev = prev;
        e5->Twin = e0;
        e5->Mark = false;

        next->Next = e3;
        next->Prev = e2;
        next->Mark = false;

        prev->Next = e5;
        prev->Prev = e4;
        prev->Mark = false;

        tri->Next = e1;
        tri->Prev = e0;
        tri->Mark = false;

        tri->recursive_swap();
        next->recursive_swap();
        prev->recursive_swap();

        return InsertPointResult::Success;
    } else {
        throw std::exception(); // TODO: No triangle could be found containing circumcenter and no boundary edge was encroached by circumcenter
    }
}

InsertPointResult Mesh::insert_circumcenter(Edge *t) {
    Point const *p = new Point(t->circumcenter());
    return insert_point(p, t);
}

InsertPointResult Mesh::insert_midpoint(Edge *e) {
    /*
        Splits edge into two edges and creates two new triangles.
    */

    Point * p = new Point;
    Points.push_back(p);

    Curve *c;
    if (e->ConstraintCurve != nullptr) { // Constrained Edge
        Vertex *v = new Vertex;
        c = e->ConstraintCurve->split(v, 0.5);

        p->X = v->x();
        p->Y = v->y();
    } else { // Unconstrained Edge
        c = nullptr;
        *(p) = Point((e->Node->X + e->Next->Node->X) / 2.0, (e->Node->Y + e->Next->Node->Y) / 2.0);
    }

    if (e == e->Twin) { // Boundary Edge
        Edge * e0 = new Edge;
        Edges.push_back(e0);

        Edge * e1 = new Edge;
        Edges.push_back(e1);

        Edge * e2 = new Edge;
        Edges.push_back(e2);

        // Save edges for swapping
        Edge * next = e->Next;
        Edge * prev = e->Prev;

        // Handle constraint curves
        e2->Orientation = e->Orientation;
        if (e->Orientation) {
            e2->ConstraintCurve = c;
        } else {
            e2->ConstraintCurve = e->ConstraintCurve;
            e->ConstraintCurve = c;
        }

        // Construct edges
        e0->Node = e->Prev->Node;
        e0->Next = e2;
        e0->Prev = e->Next;
        e0->Twin = e1;
        e0->Mark = false;

        e1->Node = p;
        e1->Next = e->Prev;
        e1->Prev = e;
        e1->Twin = e0;
        e1->Mark = false;

        e2->Node = p;
        e2->Next = e->Next;
        e2->Prev = e0;
        e2->Twin = e2;
        e2->Mark = false;

        e->Next->Next = e0;
        e->Next->Prev = e2;
        e->Next->Mark = false;

        e->Prev->Prev = e1;
        e->Prev->Mark = false;

        e->Next = e1;
        e->Mark = false;

        // Recursive swap
        e0->recursive_swap();
        next->recursive_swap();
        prev->recursive_swap();
    } else { // Interior Edge
        Edge * e0 = new Edge;
        Edges.push_back(e0);

        Edge * e1 = new Edge;
        Edges.push_back(e1);

        Edge * e2 = new Edge;
        Edges.push_back(e2);

        Edge * e3 = new Edge;
        Edges.push_back(e3);

        Edge * e4 = new Edge;
        Edges.push_back(e4);

        Edge * e5 = new Edge;
        Edges.push_back(e5);

        // Save Next/Prev for swapping
        Edge * this_next = e->Next;
        Edge * this_prev = e->Prev;
        Edge * twin_next = e->Twin->Next;
        Edge * twin_prev = e->Twin->Prev;

        // Handle constraint curves
        e1->Orientation = e->Orientation;
        e4->Orientation = !e->Orientation;
        e->Twin->Orientation = !e->Orientation;
        if (e->Orientation) {
            e1->ConstraintCurve = c;
        } else {
            e1->ConstraintCurve = e->ConstraintCurve;
            e->ConstraintCurve = c;
        }
        e->Twin->ConstraintCurve = e1->ConstraintCurve;
        e4->ConstraintCurve = e->ConstraintCurve;

        // Construct Edges
        e0->Node = p;
        e0->Next = e->Prev;
        e0->Prev = e;
        e0->Twin = e2;
        e0->Mark = false;

        e1->Node = p;
        e1->Next = e->Next;
        e1->Prev = e2;
        e1->Twin = e->Twin;
        e1->Mark = false;

        e2->Node = e->Prev->Node;
        e2->Next = e1;
        e2->Prev = e->Next;
        e2->Twin = e0;
        e2->Mark = false;

        e3->Node = p;
        e3->Next = e->Twin->Prev;
        e3->Prev = e->Twin;
        e3->Twin = e5;
        e3->Mark = false;

        e4->Node = p;
        e4->Next = e->Twin->Next;
        e4->Prev = e5;
        e4->Twin = e;
        e4->Mark = false;

        e5->Node = e->Twin->Prev->Node;
        e5->Next = e4;
        e5->Prev = e->Twin->Next;
        e5->Twin = e3;
        e5->Mark = false;

        e->Twin->Next->Next = e5;
        e->Twin->Next->Prev = e4;
        e->Twin->Next->Mark = false;

        e->Twin->Prev->Prev = e3;
        e->Twin->Prev->Mark = false;

        e->Twin->Next = e3;
        e->Twin->Twin = e1;
        e->Twin->Mark = false;

        e->Next->Next = e2;
        e->Next->Prev = e1;
        e->Twin->Mark = false;

        e->Prev->Prev = e0;
        e->Prev->Mark = false;

        e->Next = e0;
        e->Twin = e4;
        e->Mark = false;

        // Recursive swap
        e0->recursive_swap();
        e3->recursive_swap();
        this_next->recursive_swap();
        this_prev->recursive_swap();
        twin_next->recursive_swap();
        twin_prev->recursive_swap();
    }

    return InsertPointResult::Midpoint;
}

// Mesh Creation
void Mesh::create() {
    create_boundary_polygon();
    triangulate_boundary_polygon();
    get_triangles();

    insert_internal_boundaries();
    get_triangles();
}

// Algorithm Components
void Mesh::create_boundary_polygon() {
    // Create input edges
    Edges.reserve(Boundary->size());
    Points.reserve(Boundary->size());
    for (size_t i = 0; i != Boundary->size(); ++i) {
        Edges.push_back(new Edge(Boundary->curve(i)->clone(), Boundary->orientation(i))); //clone() to prevent alteration of input Contour when Edge is split
        Points.push_back(Edges[i]->Node);
    }

    // Set Next/Prev edges from Boundary
    size_t Ne = 0;
    for (size_t i = 0; i != Edges.size(); ++i) {
        size_t j = (i + 1) % Edges.size();
        Edges[i]->Next = Edges[j];
        Edges[j]->Prev = Edges[i];
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
                    if (Edges[i]->length() > Edges[j]->length()) {
                        Edges[i]->split_edge(Points, Edges);
                    } else {
                        Edges[j]->split_edge(Points, Edges);
                    }
                }
            }
        }
    }
}

void Mesh::triangulate_boundary_polygon() {
    Edges.reserve(3 * num_points());

    Edge * ep = Edges[0];
    while (ep != ep->Next->Next->Next) {
        if (ep->is_protruding()) {
            Edge * e0 = new Edge;
            Edge * e1 = new Edge;

            // Edge of new triangle
            e0->Node = ep->Next->Node;
            e0->Prev = ep;
            e0->Next = ep->Prev;
            e0->Twin = e1;

            // Twin edge, part of new polygonal boundary
            e1->Node = ep->Prev->Node;
            e1->Next = ep->Next;
            e1->Prev = ep->Prev->Prev;
            e1->Twin = e0;

            // Store new edges
            Edges.push_back(e0);
            Edges.push_back(e1);

            // Update polygonal boundary
            ep->Next->Prev = e1;
            ep->Next = e0;
            ep->Prev->Prev->Next = e1;
            ep->Prev->Prev = e0;

            // Next edge
            ep = e1->Next;
        } else {
            ep = ep->Next;
        }
    }

    // Edge swap to make triangulation Delaunay
    for (size_t i = 0; i < Edges.size(); ++i) {
        Edges[i]->recursive_swap();
    }

    split_encroached_edges();
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
        Point const *p = new Point(interior[i]->start());
        LocateTriangleResult result = locate_triangle(p);
        if (result == LocateTriangleResult::Interior) {
            while (insert_point(p) == InsertPointResult::Midpoint);
        } else {
            delete p;
        }

        // Insert end point
        p = new Point(interior[i]->end());
        result = locate_triangle(p);
        if (result == LocateTriangleResult::Interior) {
            while (insert_point(p) == InsertPointResult::Midpoint);
        } else {
            delete p;
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

            LocateTriangleResult result = locate_triangle(&p0, e);

            if (result != LocateTriangleResult::Point) {
                throw std::exception();
            }

            if (e->is_attached(&p1, e)) {
                e->ConstraintCurve = queue.back();
                e->Orientation = true;

                e->Twin->ConstraintCurve = queue.back();
                e->Twin->Orientation = false;

                queue.pop_back();
            } else {
                Vertex *v = new Vertex;
                queue.push_back(queue.back()->split(v, 0.5));

                Point const *p = new Point(queue.back()->start());

                while (insert_point(p) == InsertPointResult::Midpoint);
            }
        }
    }

    split_encroached_edges();
}

void Mesh::split_encroached_edges() {
    bool any_split = true;
    while (any_split) {
        any_split = false;
        for (size_t i = 0; i != Edges.size(); ++i) {
            if (Edges[i]->constraint_curve() != nullptr) {
                if (Edges[i]->is_encroached(Edges[i]->prev()->base())) {
                    any_split = true;
                    insert_midpoint(Edges[i]);
                }
            }
        }
    }
}

// Misc
void Mesh::mark_triangles() {
    for (size_t i = 0; i < Edges.size(); ++i) {
        Edges[i]->Mark = true;
    }

    for (size_t i = 0; i < Edges.size(); ++i) {
        Edges[i]->recursive_mark();
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

// Save
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
        Point const *v0 = Edges[i]->Node;
        Point const *v1 = Edges[i]->Next->Node;
        Point const *v2 = Edges[i]->Next->Next->Node;
        fs << v0->X << ',' << v1->X << ',' << v2->X << ',' << v0->Y << ',' << v1->Y << ',' << v2->Y << '\n';
    }

    fs.close();
}

// Refinement
bool Mesh::refine() {
    std::vector<double> radii;
    std::vector<double> quality;
    std::vector<size_t> index;

    // #TODO, Loop until quality is satisfied
    element_quality(Triangles, radii, quality);
    sort_permutation(quality, index);
    size_t N = Triangles.size();

    refine_once(index, radii, quality);
    size_t M = Triangles.size();

    while (M > N) {
        N = M;
        element_quality(Triangles, radii, quality);
        sort_permutation(quality, index);
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
    sort_permutation(quality, index);
    refine_once(index, radii, quality);

    return edges_are_valid();
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

bool Mesh::edges_are_valid() {
    bool result = true;
    Edge const *e;

    for (size_t i = 0; i < size_edges(); ++i) {
        e = edge(i);

        if (e != e->next()->prev()) {
            result = false;
            break;
        }
        if (e != e->prev()->next()) {
            result = false;
            break;
        }
        if (e != e->twin()->twin()) {
            result = false;
            break;
        }

        if (!(e->twin() == e)) {
            if (e->node() != e->twin()->next()->node()) {
                result = false;
                break;
            }
            if (e->constraint_curve() != e->twin()->constraint_curve()) {
                result = false;
                break;
            }
            if (e->constraint_curve() != nullptr) {
                if (e->orientation() == e->twin()->orientation()) {
                    result = false;
                    break;
                }
            }

            if (e->node() == e->twin()->node()) {
                result = false;
                break;
            }
        }

        if (e->constraint_curve() != nullptr) {
            if (e->orientation()) {
                if (*e->base() != *e->constraint_curve()->start()) {
                    result = false;
                    break;
                }
                if (*e->tip() != *e->constraint_curve()->end()) {
                    result = false;
                    break;
                }
            } else {
                if (*e->base() != *e->constraint_curve()->end()) {
                    result = false;
                    break;
                }
                if (*e->tip() != *e->constraint_curve()->start()) {
                    result = false;
                    break;
                }
            }
        }
    }

    return result;
}