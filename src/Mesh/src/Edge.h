#ifndef OERSTED_EDGE_H
#define OERSTED_EDGE_H

#include "Mesh.h"

class Edge {
public:
    friend class Mesh;

    friend bool are_intersecting(const Edge *e0, const Edge *e1);

    // Constructors
    Edge() : Node(nullptr), Next(nullptr), Prev(nullptr), Twin(this), ConstraintCurve(nullptr), Orientation(true), Mark(false) {};

    Edge(Point &v, Edge &n, Edge &p, Edge &t) : Node(&v), Next(&n), Prev(&p), Twin(&t), ConstraintCurve(nullptr), Orientation(true), Mark(false) {};

    Edge(Curve *c, bool Orientation);

    Edge(Curve *c, bool Orientation, const Point *v);

    // Accessors
    const Point *node() const { return Node; };

    const Point *base() const { return Node; };

    const Point *tip() const { return (next() == nullptr ? twin()->base() : next()->base()); };

    const Edge *next() const { return Next; };

    const Edge *twin() const { return Twin; };

    const Edge *prev() const { return Prev; };

    const Curve *constraint_curve() const { return ConstraintCurve; };

    bool orientation() const { return Orientation; };

    bool mark() const { return Mark; };

    bool operator==(const Edge &e) const {
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

    bool is_encroached(const Point *p) const;

    bool is_attached(const Point *p, Edge *&e) const;

    bool swap();

    bool recursive_swap();

    void split_edge(std::vector<const Point *> &verts, std::vector<Edge *> &edges);

    void recursive_mark();

protected:
    const Point *Node;      //Start of edge
    Edge *Next;             //In triangle
    Edge *Twin;             //Adjacent triangle
    Edge *Prev;             //In triangle

    Curve *ConstraintCurve; //==nullptr if unconstrained
    bool Orientation;       //undefined if unconstrained

    bool Mark;              // Auxillary variable for mesh refinement
};

#endif //OERSTED_EDGE_H
