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
    friend class Edge;

public:
    double MinimumElementQuality = 0.0;
    double MinimumElementSize = 0.0;
    double MaximumElementSize = DBL_MAX;

    Mesh(Sketch &s);

    bool are_intersecting(Edge const *e0, Edge const *e1) const;

    bool edges_are_valid() const;

    bool in_triangle(Point const p, Edge const *e) const;

    bool is_encroached(Edge const *e, Point const p) const;

    bool is_locally_optimal(Edge const *e) const;

    bool is_protruding(Edge const *e) const;

    bool is_valid(Edge const *e) const;

    bool refine();

    bool refine_once();

    double circumradius(Edge const *e) const;

    double length(Edge const *e) const;

    double shortest_edge_length(Edge const *e) const;

    size_t size_points() const { return Points.size(); };

    size_t size_edges() const { return Edges.size(); };

    size_t size_triangles() const { return Triangles.size(); };

    size_t num_points() const { return Points.size(); };

    size_t num_edges() const;

    size_t num_triangles() const { return Triangles.size(); };

    void create();

    void delete_me(); // TODO: refactor to non-pointer version

    void save_as(std::string path, std::string file_name) const;

    Point circumcenter(Edge const *e) const;

    Point const point(size_t i) const { return Points[i]; };

    Edge const *edge(size_t i) const { return Edges[i]; };

    Edge const *triangle(size_t i) const { return Triangles[i]; };

    LocateTriangleResult locate_triangle(Point const p, Edge const *&e) const;

    LocateTriangleResult locate_triangle(Point const p) const {
        Edge const *e = Edges.back();
        return locate_triangle(p, e);
    };

    InsertPointResult insert_point(Point const p) { return insert_point(p, Edges.back()); };

protected:
    Contour const *Boundary;
    std::vector<Curve const *> Curves;
    std::vector<Contour const *> Contours;
    std::vector<Point> Points;
    std::vector<Edge *> Edges;
    std::vector<Edge *> Triangles;

private:
    bool find_attached(Edge *&e_out, Point const p) const;

    bool recursive_swap(Edge *e) const;

    bool swap(Edge *&e0) const;

    void add_edge(Edge *&e) {
        e->Self = Edges.size();
        Edges.push_back(e);
    };

    void create_boundary_polygon();

    void element_quality(std::vector<Edge *> &triangle, std::vector<double> &radii, std::vector<double> &quality);

    void get_triangles();

    void insert_internal_boundaries();

    void mark_triangles();

    void recursive_mark(Edge *e) const;

    void refine_once(std::vector<size_t> index, std::vector<double> circumradius, std::vector<double> quality);

    void sort_permutation_ascending(std::vector<double> &value, std::vector<size_t> &index) const;

    void sort_permutation_descending(std::vector<double> &values, std::vector<size_t> &index) const;

    void split_edge(Edge *e);

    void split_encroached_edges();

    void triangulate_boundary_polygon();

    LocateTriangleResult locate_triangle(Point const p, Edge *&e) const;

    InsertPointResult insert_circumcenter(Edge *e);

    InsertPointResult insert_point(Point const p, Edge *e);

    InsertPointResult insert_midpoint(Edge *e);
};

#endif //OERSTED_MESH_H
