#ifndef OERSTED_MESH_H
#define OERSTED_MESH_H

#include "Sketch.hpp"
#include "Point.h"
#include "Edge.h"

#include <fstream>
#include <algorithm>
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

    void delete_me(); // TODO: refactor to non-pointer version

    void create();

    Point circumcenter(Edge const *e) const;

    Point const point(size_t i) const { return Points[i]; };

    Edge const *edge(size_t i) const { return Edges[i]; };

    Edge const *triangle(size_t i) const { return Triangles[i]; };

    size_t size_points() const { return Points.size(); };

    size_t size_edges() const { return Edges.size(); };

    size_t size_triangles() const { return Triangles.size(); };

    size_t num_points() const { return Points.size(); };

    size_t num_edges() const;

    size_t num_triangles() const { return Triangles.size(); };

    void save_as(std::string path, std::string file_name) const;

    LocateTriangleResult locate_triangle(Point const &p, Edge *&e) const;

    LocateTriangleResult locate_triangle(Point const &p, Edge const *&e) const;

    LocateTriangleResult locate_triangle(Point const p) const {
        Edge const *e = Edges.back();
        return locate_triangle(p, e);
    };

    InsertPointResult insert_point(Point const &p, Edge *e);

    InsertPointResult insert_point(Point const p) { return insert_point(p, Edges.back()); };

    InsertPointResult insert_circumcenter(Edge *e);

    InsertPointResult insert_midpoint(Edge *e);

    bool refine();

    bool refine_once();

    void refine_once(std::vector<size_t> index, std::vector<double> circumradius, std::vector<double> quality);

    bool edges_are_valid() const;

    double circumradius(Edge const *e) const;

    double length(Edge const *e) const;

    double shortest_edge_length(Edge const *e) const;

    bool is_protruding(Edge const *e) const;

    bool recursive_swap(Edge *edge) const;

    bool is_locally_optimal(Edge const *edge) const;

    void split_edge(Edge *edge);

    bool is_encroached(Edge const *e, Point const p2) const;

    bool is_attached(Edge *&edge_out, Point const &p) const;

protected:
    Contour const *Boundary;
    std::vector<Curve const *> Curves;
    std::vector<Contour const *> Contours;
    std::vector<Point> Points;
    std::vector<Edge *> Edges;
    std::vector<Edge *> Triangles;

private:
    void create_boundary_polygon();

    void triangulate_boundary_polygon();

    void insert_internal_boundaries();

    void split_encroached_edges();

    // Misc
    void mark_triangles();

    void get_triangles();
};

#endif //OERSTED_MESH_H
