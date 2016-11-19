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

class DartConstraint {
public:
    DartConstraint() : S0(DBL_MAX), S1(DBL_MAX), ConstraintCurve(nullptr) {};

    DartConstraint(double s0, double s1, Curve const *cc) : S0(s0), S1(s1), ConstraintCurve(cc) {};

    double S0;
    double S1;
    Curve const *ConstraintCurve;
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

    bool is_constrained(size_t ei) const { return Edges[ei].Constraint != 0; };

    bool is_encroached(Point const p, size_t ei) const;

    bool is_locally_optimal(size_t ei) const;

    bool is_protruding(size_t ei) const;

    bool is_valid(size_t ei) const;

    bool orientation(size_t ei) const { return Edges[ei].Orientation; };

    bool refine();

    bool refine_once();

    double circumradius(size_t ei) const;

    double length(size_t ei) const;

    double shortest_edge_length(size_t ei) const;

    size_t next(size_t ei) const { return Edges[ei].Next; };

    size_t node(size_t ei) const { return Edges[ei].Node; };

    size_t node(Edge const e) const { return e.Node; };

    size_t num_points() const { return Points.size(); };

    size_t num_edges() const;

    size_t num_triangles() const { return Triangles.size(); };

    size_t prev(size_t ei) const { return Edges[ei].Prev; };

    size_t size_points() const { return Points.size(); };

    size_t size_edges() const { return Edges.size(); };

    size_t size_triangles() const { return Triangles.size(); };

    size_t twin(size_t ei) const { return Edges[ei].Twin; };

    void create();

    void save_as(std::string path, std::string file_name) const;

    DartConstraint const constraint(size_t ei) const { return Constraints[Edges[ei].Constraint]; };

    Curve const *constraint_curve(size_t ei) const { return Constraints[Edges[ei].Constraint].ConstraintCurve; };

    Point circumcenter(size_t ei) const;

    Point const base(Edge const e) const { return Points[e.Node]; };

    Point const base(size_t ei) const { return Points[node(ei)]; };

    Point const point(size_t i) const { return Points[i]; };

    Point const point(Edge const e) const { return Points[e.Node]; };

    Point const tip(Edge const e) const { return Points[next(e).Node]; };

    Point const tip(size_t ei) const { return Points[node(next(ei))]; };

    Edge const edge(size_t i) const { return Edges[i]; };

    Edge const next(Edge const e) const { return Edges[e.Next]; };

    Edge const prev(Edge const e) const { return Edges[e.Prev]; };

    Edge const twin(Edge const e) const { return Edges[e.Twin]; };

    Edge const triangle(size_t i) const { return Edges[Triangles[i]]; };

    LocateTriangleResult locate_triangle(Point const p, size_t &ei) const;

    LocateTriangleResult locate_triangle(Point const p) const {
        size_t ei = Edges.size() - 1;
        return locate_triangle(p, ei);
    };

    InsertPointResult insert_point(Point const p) { return insert_point(p, Edges.size() - 1); };

protected:
    std::shared_ptr<Contour> Boundary;
    std::vector<Curve const *> Curves;
    std::vector<std::shared_ptr<Contour>> Contours;

    std::vector<Point> Points;
    std::vector<Edge> Edges;
    std::vector<DartConstraint> Constraints;

    std::vector<size_t> Triangles;

private:
    bool find_attached(Point const p, size_t &ei);

    bool recursive_swap(size_t ei);

    bool swap(size_t ei);

    size_t new_edges(size_t num_new) {
        for (size_t i = 0; i != num_new; ++i) {
            Edges.push_back(Edge(Edges.size()));
            Edges.back().Constraint = 0;
        }
        return Edges.size();
    }

    Edge &new_edge(size_t p, size_t c, bool dir) {
        Edges.push_back(Edge(p, Edges.size(), c, dir));
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
