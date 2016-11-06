#ifndef OERSTED_EDGE_H
#define OERSTED_EDGE_H

#include "Point.h"

class Mesh;

class Edge {
public:
    friend class Mesh;

    friend bool are_intersecting(Edge const *e0, Edge const *e1, Mesh const &m);

    // Constructors
    Edge() : Node(SIZE_MAX), Next(nullptr), Prev(nullptr), Twin(this), ConstraintCurve(nullptr), Orientation(true), Mark(false) {};

    Edge(size_t v, Edge &n, Edge &p, Edge &t) : Node(v), Next(&n), Prev(&p), Twin(&t), ConstraintCurve(nullptr), Orientation(true), Mark(false) {};

    Edge(Curve *c, bool Orientation, size_t v);

    // Accessors
    size_t node() const { return Node; };

    size_t base() const { return Node; };

    size_t tip() const { return (next() == nullptr ? twin()->base() : next()->base()); };

    Edge const *next() const { return Next; };

    Edge const *twin() const { return Twin; };

    Edge const *prev() const { return Prev; };

    Curve const *constraint_curve() const { return ConstraintCurve; };

    bool orientation() const { return Orientation; };

    bool mark() const { return Mark; };

    bool operator==(Edge const &e) const {
        return (node() == e.node()) && (constraint_curve() == e.constraint_curve()) && (twin()->node() == e.twin()->node());
    };

    bool is_constrained() const { return (ConstraintCurve != nullptr); };

    bool is_valid() const;

    bool swap();

    void recursive_mark();

protected:
    size_t Node;            //Start of edge
    Edge *Next;             //In triangle
    Edge *Twin;             //Adjacent triangle
    Edge *Prev;             //In triangle

    Curve *ConstraintCurve; //==nullptr if unconstrained
    bool Orientation;       //undefined if unconstrained

    bool Mark;              // Auxillary variable for mesh refinement
};

#endif //OERSTED_EDGE_H
