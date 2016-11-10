#ifndef OERSTED_EDGE_H
#define OERSTED_EDGE_H

#include "Point.h"

class Mesh;

class Edge {
public:
    friend class Mesh;

    friend bool are_intersecting(Edge const *e0, Edge const *e1, Mesh const &m);

    // Constructors
    Edge() : Node(SIZE_MAX), Self(SIZE_MAX), Next(SIZE_MAX), Prev(SIZE_MAX), Twin(SIZE_MAX), ConstraintCurve(nullptr), Orientation(true), Mark(false) {};

    Edge(size_t v, size_t s, Edge &n, Edge &p, Edge &tw) : Node(v), Self(s), Next(n.Self), Prev(p.Self), Twin(tw.Self), ConstraintCurve(nullptr), Orientation(true), Mark(false) {};

    Edge(Curve *c, bool Orientation, size_t v, size_t s);

    // Accessors
    size_t node() const { return Node; };

    size_t base() const { return Node; };

    size_t tip(Mesh const &mesh) const;

    size_t self() const { return Self; };

    size_t next() const { return Next; };

    size_t twin() const { return Twin; };

    size_t prev() const { return Prev; };

    Curve const *constraint_curve() const { return ConstraintCurve; };

    bool orientation() const { return Orientation; };

    bool mark() const { return Mark; };

    bool operator==(Edge const &e) const {
        return (Node == e.Node &&
                Self == e.Self &&
                Next == e.Next &&
                Twin == e.Twin &&
                Prev == e.Prev &&
                ConstraintCurve == e.ConstraintCurve &&
                Orientation == e.Orientation);
    };

    bool is_constrained() const { return (ConstraintCurve != nullptr); };

protected:
    size_t Node;            //Start of edge
    size_t Self;            //This triangle
    size_t Next;            //In triangle
    size_t Twin;            //Adjacent triangle
    size_t Prev;            //In triangle

    Curve *ConstraintCurve; //==nullptr if unconstrained
    bool Orientation;       //undefined if unconstrained

    bool Mark;              // Auxillary variable for mesh refinement
};

#endif //OERSTED_EDGE_H
