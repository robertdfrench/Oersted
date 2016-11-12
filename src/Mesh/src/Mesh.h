#ifndef OERSTED_MESH_H
#define OERSTED_MESH_H

#include "Sketch.hpp"
#include "Point.h"
#include "Edge.h"

#include <algorithm>
#include <fstream>
#include <numeric>

enum class LocateTriangleResult {
    Interior, Exterior, Boundary, Point // TODO: Enumerate cases in function when triangle is located by a point near the boundary
};

enum class InsertPointResult {
    Success, Midpoint, Duplicate, Failure
};

class Mesh {
public:
    double MinimumElementQuality = 0.0;
    double MinimumElementSize = 0.0;
    double MaximumElementSize = DBL_MAX;

    Mesh(Sketch &s);

    bool are_intersecting(size_t ei, size_t ej) const;

    bool edges_are_valid() const;

    bool in_triangle(Point const p, size_t ei) const;

    bool is_constrained(size_t ei) const { return Edges[ei]->is_constrained(); };

    bool is_encroached(Point const p, size_t ei) const;

    bool is_locally_optimal(size_t ei) const;

    bool is_protruding(size_t ei) const;

    bool is_valid(size_t ei) const;

    bool orientation(size_t ei) const { return Edges[ei]->Orientation; };

    bool refine();

    bool refine_once();

    double circumradius(size_t ei) const;

    double length(size_t ei) const;

    double shortest_edge_length(size_t ei) const;

    size_t next(size_t ei) const { return Edges[ei]->Next; };

    size_t node(size_t ei) const { return Edges[ei]->Node; };

    size_t node(Edge const *e) const { return e->Node; };

    size_t num_points() const { return Points.size(); };

    size_t num_edges() const;

    size_t num_triangles() const { return Triangles.size(); };

    size_t prev(size_t ei) const { return Edges[ei]->Prev; };

    size_t size_points() const { return Points.size(); };

    size_t size_edges() const { return Edges.size(); };

    size_t size_triangles() const { return Triangles.size(); };

    size_t twin(size_t ei) const { return Edges[ei]->Twin; };

    void create();

    void delete_me(); // TODO: refactor to non-pointer version

    void save_as(std::string path, std::string file_name) const;

    Curve const *constraint_curve(size_t ei) const { return Edges[ei]->ConstraintCurve; };

    Point circumcenter(size_t ei) const;

    Point const base(Edge const *e) const { return Points[e->Node]; };

    Point const base(size_t ei) const { return Points[node(ei)]; };

    Point const point(size_t i) const { return Points[i]; };

    Point const point(Edge const *e) const { return Points[e->Node]; };

    Point const tip(Edge const *e) const { return Points[next(e)->Node]; };

    Point const tip(size_t ei) const { return Points[node(next(ei))]; };

    Edge const *edge(size_t i) const { return Edges[i]; };

    Edge const *next(Edge const *e) const { return Edges[e->Next]; };

    Edge *&next(Edge *e) { return Edges[e->Next]; };

    Edge const *prev(Edge const *e) const { return Edges[e->Prev]; };

    Edge *&prev(Edge *e) { return Edges[e->Prev]; };

    Edge const *twin(Edge const *e) const { return Edges[e->Twin]; };

    Edge *&twin(Edge *e) { return Edges[e->Twin]; };

    Edge const *triangle(size_t i) const { return Edges[Triangles[i]]; };

    LocateTriangleResult locate_triangle(Point const p, size_t &ei) const;

    LocateTriangleResult locate_triangle(Point const p) const {
        //Edge const *e = Edges.back();
        size_t ei = Edges.size() - 1;
        return locate_triangle(p, ei);
    };

    InsertPointResult insert_point(Point const p) { return insert_point(p, Edges.size() - 1); };

protected:
    Contour const *Boundary;
    std::vector<Curve const *> Curves;
    std::vector<Contour const *> Contours;
    std::vector<Point> Points;
    std::vector<Edge *> Edges;
    std::vector<size_t> Triangles;

private:
    bool find_attached(Point const p, size_t &ei);

    bool recursive_swap(size_t ei);

    bool swap(size_t ei);

    Edge *&new_edge() {
        Edges.push_back(new Edge(Edges.size()));
        return Edges.back();
    }

    Edge *&new_edge(size_t p, Curve *c, bool dir) {
        Edges.push_back(new Edge(p, Edges.size(), c, dir));
        return Edges.back();
    }

    void create_boundary_polygon();

    void element_quality(std::vector<double> &radii, std::vector<double> &quality);

    void get_triangles();

    void insert_internal_boundaries();

    void mark_triangles();

    void refine_once(std::vector<size_t> index, std::vector<double> circumradius, std::vector<double> quality);

    void sort_permutation_ascending(std::vector<double> &value, std::vector<size_t> &index) const;

    void sort_permutation_descending(std::vector<double> &values, std::vector<size_t> &index) const;

    void split_edge(size_t ei);

    void split_encroached_edges();

    void triangulate_boundary_polygon();

    InsertPointResult insert_circumcenter(size_t ei);

    InsertPointResult insert_point(Point const p, size_t ei);

    InsertPointResult insert_midpoint(size_t ei);
};

#endif //OERSTED_MESH_H
