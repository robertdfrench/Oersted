#ifndef OERSTED_EDGE_H
#define OERSTED_EDGE_H

#include "Mesh.h"

class Edge {
public:
    friend class Mesh;

    friend bool are_intersecting(Edge const *e0, Edge const *e1);

    // Constructors
    Edge() : Node(nullptr), Next(nullptr), Prev(nullptr), Twin(this), ConstraintCurve(nullptr), Orientation(true), Mark(false) {};

    Edge(Point &v, Edge &n, Edge &p, Edge &t) : Node(&v), Next(&n), Prev(&p), Twin(&t), ConstraintCurve(nullptr), Orientation(true), Mark(false) {};

    Edge(Curve *c, bool Orientation);

    Edge(Curve *c, bool Orientation, Point const *v);

    // Accessors
    Point const *node() const { return Node; };

    Point const *base() const { return Node; };

    Point const *tip() const { return (next() == nullptr ? twin()->base() : next()->base()); };

    Edge const *next() const { return Next; };

    Edge const *twin() const { return Twin; };

    Edge const *prev() const { return Prev; };

    Curve const *constraint_curve() const { return ConstraintCurve; };

    bool orientation() const { return Orientation; };

    bool mark() const { return Mark; };

    bool operator==(Edge const &e) const {
        return (node() == e.node()) && (constraint_curve() == e.constraint_curve()) && (twin()->node() == e.twin()->node());
    };

    Point circumcenter() const;

    double circumradius() const;

    double length() const;

    double shortest_edge_length() const;

    bool is_protruding() const;

    bool is_constrained() const { return (ConstraintCurve != nullptr); };

    bool is_locally_optimal() const;

    bool is_valid() const;

    bool is_encroached(Point const *p) const;

    bool is_attached(Point const *p, Edge *&e) const;

    bool swap();

    bool recursive_swap();

    void split_edge(std::vector<Point const *> &verts, std::vector<Edge *> &edges);

    void recursive_mark();

protected:
    Point const *Node;      //Start of edge
    Edge *Next;             //In triangle
    Edge *Twin;             //Adjacent triangle
    Edge *Prev;             //In triangle

    Curve *ConstraintCurve; //==nullptr if unconstrained
    bool Orientation;       //undefined if unconstrained

    bool Mark;              // Auxillary variable for mesh refinement
};

#endif //OERSTED_EDGE_H
