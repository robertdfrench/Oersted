#ifndef OERSTED_MESH_H
#define OERSTED_MESH_H

#include "Sketch.hpp"

#include <fstream>
#include <algorithm>
#include <numeric>

enum class LocateTriangleResult {
    Interior, Exterior, Edge, Point
};

enum class InsertPointResult {
    Success, Midpoint, Duplicate, Failure
};

class Point { // #TODO, Replace Verticies member in Mesh with lighterweight Node class
public:

    Point() : X{0.0}, Y{0.0} {};

    Point(double x, double y) : X{x}, Y{y} {};

    Point(const Vertex &v) : X{v.x()}, Y{v.y()} {};

    Point(const Vertex *v) : X{v->x()}, Y{v->y()} {};

    //double W; // Nurbs weight?
    double X;
    double Y;
    //double Z; // 3-Dimensions?

    //Vertex* V = nullptr; Keep track of constrained nodes?

    bool operator==(const Point &p) const { return (X == p.X) && (Y == p.Y); };

    bool operator==(const Vertex &v) const { return (X == v.x()) && (Y == v.y()); };

    bool operator!=(const Point &p) const { return (X != p.X) || (Y != p.Y); };

    bool operator!=(const Vertex &v) const { return (X != v.x()) && (Y != v.y()); };
};

class Edge;

class Mesh {
public:
    // Global Mesh Quality Constraints
    double MinimumElementQuality = 0.0;
    double MinimumElementSize = 0.0;
    double MaximumElementSize = DBL_MAX;

    // Constructors
    Mesh(Sketch &s);

    void create();

    // Accessors
    const Point *point(size_t i) const { return Points[i]; };

    const Edge *edge(size_t i) const { return Edges[i]; };

    const Edge *triangle(size_t i) const { return Triangles[i]; };

    size_t size_points() const { return Points.size(); };

    size_t size_edges() const { return Edges.size(); };

    size_t size_triangles() const { return Triangles.size(); };

    size_t num_points() const { return Points.size(); };

    size_t num_edges() const;

    size_t num_triangles() const { return Triangles.size(); };

    // Save
    void save_as(std::string path, std::string file_name) const;

    // Topological Queries
    LocateTriangleResult locate_triangle(const Point *p, Edge *&e) const;

    LocateTriangleResult locate_triangle(const Point *p, const Edge *&e) const;

    LocateTriangleResult locate_triangle(const Point *p) const {
        const Edge *e = Edges.back();
        return locate_triangle(p, e);
    };

    // Point Insertion
    InsertPointResult insert_point(const Point *p, Edge *e);

    InsertPointResult insert_point(const Point *p) { return insert_point(p, Edges.back()); };

    InsertPointResult insert_circumcenter(Edge *e);

    InsertPointResult insert_midpoint(Edge *e);

    // Refinement
    void refine();

    void refine_once();

    void refine_once(std::vector<size_t> index, std::vector<double> circumradius, std::vector<double> quality);

protected:
    const Contour *Boundary;
    std::vector<const Curve *> Curves;
    std::vector<const Contour *> Contours;
    std::vector<const Point *> Points;
    std::vector<Edge *> Edges;
    std::vector<Edge *> Triangles;

private:
    // Algorithm Components
    void create_boundary_polygon();

    void triangulate_boundary_polygon();

    void insert_internal_boundaries();

    void split_encroached_edges();

    // Misc
    void mark_triangles();

    void get_triangles();
};

#endif //OERSTED_MESH_H
