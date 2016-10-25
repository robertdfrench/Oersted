#include "Mesh.hpp"

Mesh::Mesh(Sketch &s) {
    Boundary = s.boundary();

    for (size_t i = 0; i != s.size_curves(); ++i) {
        Curves.push_back(s.curve(i));
    }

    for (size_t i = 0; i != s.size_contours(); ++i) {
        Contours.push_back(s.contour(i));
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

// Topological Queries
LocateTriangleResult Mesh::locate_triangle(const Point *p, Edge *&e) const {
    const Edge *ec = e;

    LocateTriangleResult result = locate_triangle(p, ec);

    e = const_cast<Edge *>(ec);

    return result;
}

LocateTriangleResult Mesh::locate_triangle(const Point *p, const Edge *&e) const {
    const Edge *start = e;
    double xp = p->X;
    double yp = p->Y;

    while (true) {
        // Verticies
        const Point *p0 = e->base();
        const Point *p1 = e->tip();

        if (*p == *p0) {
            return LocateTriangleResult::Point;
        } else if (*p == *p1) {
            e = e->next();
            return LocateTriangleResult::Point;
        }

        // Signed area of triangle
        double dx0p = p0->X - xp;
        double dy0p = p0->Y - yp;
        double dx1p = p1->X - xp;
        double dy1p = p1->Y - yp;
        double dx01 = p0->X - p1->X;
        double dy01 = p0->Y - p1->Y;

        double area = dx0p * dy1p - dx1p * dy0p;
        double tol = FLT_EPSILON * (dx01 * dx01 + dy01 * dy01);

        if (area < -tol) {
            if (e != e->twin()) {
                start = e->twin();
                e = start->next();
            } else {
                return LocateTriangleResult::Exterior;
            }
        } else if (area < tol) {
            double dp0 = sqrt(dx0p * dx0p + dy0p * dy0p);
            double dp1 = sqrt(dx1p * dx1p + dy1p * dy1p);
            double d01 = sqrt(dx01 * dx01 + dy01 * dy01);
            double tol = FLT_EPSILON * d01;

            if (abs(dp0 + dp1 - d01) < tol) {
                return LocateTriangleResult::Edge;
            } else if (e->next() != start) {
                e = e->next();
            } else {
                throw (std::exception()); //TODO: This branch should never be reached. Supply return value"
            }
        } else if (e->next() != start) {
            e = e->next();
        } else {
            return LocateTriangleResult::Interior;
        }
    }
}

// Point Insertion
InsertPointResult Mesh::insert_point(const Point *vc, Edge *tri) {
    // Find triangle containing point
    LocateTriangleResult result = locate_triangle(vc, tri);

    // Do not insert point if it is duplicate
    if (result == LocateTriangleResult::Point) {
        return InsertPointResult::Duplicate;
    }

    // Test edges in current and adjacent triangles for encroachment
    // These are the only possible edges that are encroached due to empy circumcircle property
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

        const Point *vt = tri->Node;
        const Point *vn = next->Node;
        const Point *vp = prev->Node;

        Edge * e0 = new Edge;
        Edge * e1 = new Edge;
        Edge * e2 = new Edge;
        Edge * e3 = new Edge;
        Edge * e4 = new Edge;
        Edge * e5 = new Edge;

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
    } else if (result == LocateTriangleResult::Edge) {
        // Save edges for swapping
        Edge * next = tri->Next;
        Edge * prev = tri->Prev;
        Edge * twin = tri->Twin;
        Edge * twin_next = tri->Twin->Next;
        Edge * twin_prev = tri->Twin->Prev;

        // Get vertex pointers
        const Point *vt = tri->Node;
        const Point *vn = next->Node;
        const Point *vp = prev->Node;
        const Point *vtp = twin_prev->Node;

        // Allocate new edges
        Edge * e0 = new Edge;
        Edge * e1 = new Edge;
        Edge * e2 = new Edge;
        Edge * e3 = new Edge;
        Edge * e4 = new Edge;
        Edge * e5 = new Edge;

        // Save
        Points.push_back(vc);

        Edges.push_back(e0);
        Edges.push_back(e1);
        Edges.push_back(e2);
        Edges.push_back(e3);
        Edges.push_back(e4);
        Edges.push_back(e5);

        // Create new triangles
        e0->Node = vc;
        e0->Next = prev;
        e0->Prev = tri;
        e0->Twin = e1;
        e0->Mark = false;

        e1->Node = vp;
        e1->Next = e2;
        e1->Prev = next;
        e1->Twin = e0;
        e1->Mark = false;

        e2->Node = vc;
        e2->Next = next;
        e2->Prev = e1;
        e2->Twin = twin;
        e2->Mark = false;

        e3->Node = vc;
        e3->Next = twin_prev;
        e3->Prev = twin;
        e3->Twin = e4;
        e3->Mark = false;

        e4->Node = vtp;
        e4->Next = e5;
        e4->Prev = twin_next;
        e4->Twin = e3;
        e4->Mark = false;

        e5->Node = vc;
        e5->Next = twin_next;
        e5->Prev = e4;
        e5->Twin = tri;
        e5->Mark = false;

        twin_next->Next = e4;
        twin_next->Prev = e5;
        twin_next->Mark = false;

        twin_prev->Prev = e3;
        twin_prev->Mark = false;

        twin->Next = e3;
        twin->Twin = e2;
        twin->Mark = false;

        next->Next = e1;
        next->Prev = e2;
        next->Mark = false;

        prev->Prev = e0;
        prev->Mark = false;

        tri->Next = e0;
        tri->Twin = e5;
        tri->Mark = false;

        next->recursive_swap();
        prev->recursive_swap();
        twin_next->recursive_swap();
        twin_prev->recursive_swap();

        return InsertPointResult::Success;
    } else {
        throw (std::exception()); // TODO: No triangle could be found containing circumcenter and no edge was encroached by circumcenter
    }
}

InsertPointResult Mesh::insert_circumcenter(Edge *t) {
    const Point *p = new Point(t->circumcenter());
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
        Edges.push_back(new Edge(Boundary->curve(i)->clone(), Boundary->orientation(
                i))); //clone() to prevent alteration of input Contour when Edge is split
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

    // Create initial polygon
    size_t Ne = 0;
    for (size_t i = 0; i < Contours.size(); ++i) {
        Edge * ep = Edges[Ne];
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
        Ne += Contours[i]->size();
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
        const Point *p = new Point(interior[i]->start());
        LocateTriangleResult result = locate_triangle(p);
        if (result == LocateTriangleResult::Interior || result == LocateTriangleResult::Edge) {
            while (insert_point(p) == InsertPointResult::Midpoint);
        } else {
            delete p;
        }

        // Insert end point
        p = new Point(interior[i]->end());
        result = locate_triangle(p);
        if (result == LocateTriangleResult::Interior || result == LocateTriangleResult::Edge) {
            while (insert_point(p) == InsertPointResult::Midpoint);
        } else {
            delete p;
        }
    }

    // Insert interior curve midpoints until constraints are naturally satisified
    std::vector<Curve *> queue;
    for (size_t i = 0; i != interior.size(); ++i) {
        Edge * e = Edges.back();
        queue.push_back(interior[i]);
        while (queue.size() != 0) {
            Point p0 = queue.back()->start();
            Point p1 = queue.back()->end();

            LocateTriangleResult result = locate_triangle(&p0, e);

            if (result != LocateTriangleResult::Point) {
                throw (std::exception()); //TODO: Could not locate point within triangulation while inserting internal constraints
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

                const Point *p = new Point(queue.back()->start());

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
    fs.open(path + file_name, std::fstream::out);

    for (size_t i = 0; i < Edges.size(); i++) {
        const Point *v0 = Edges[i]->Node;
        const Point *v1 = Edges[i]->Next->Node;
        const Point *v2 = Edges[i]->Next->Next->Node;
        fs << v0->X << "," << v1->X << "," << v2->X << "," << v0->Y << "," << v1->Y << "," << v2->Y << "\n";
    }

    fs.close();
}

// Refinement
void Mesh::refine() {
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
}

void Mesh::refine_once() {
    std::vector<double> radii;
    std::vector<double> quality;
    std::vector<size_t> index;

    element_quality(Triangles, radii, quality);
    sort_permutation(quality, index);
    refine_once(index, radii, quality);
}

void Mesh::refine_once(std::vector<size_t> index, std::vector<double> radii, std::vector<double> quality) {
    for (size_t i = 0; i < Triangles.size(); ++i) {
        size_t j = index[i];
        if ((Triangles[j]->Mark) && ((radii[j] > MaximumElementSize) || (radii[j] > MinimumElementSize && quality[j] < MinimumElementQuality))) {
            //Triangles[j]->insert_circumcenter(Points, Edges);
            insert_circumcenter(Triangles[j]);
        }
    }
    get_triangles();
}