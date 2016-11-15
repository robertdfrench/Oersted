#ifndef OERSTED_EDGE_H
#define OERSTED_EDGE_H

#include "Point.h"

class Mesh;

class Edge {
public:
    friend class Mesh;

    Edge() : Node(SIZE_MAX), Self(SIZE_MAX), Next(SIZE_MAX), Twin(SIZE_MAX), Prev(SIZE_MAX), Constraint(SIZE_MAX), Orientation(true), Mark(false) {};

    Edge(size_t s) : Node(SIZE_MAX), Self(s), Next(SIZE_MAX), Twin(SIZE_MAX), Prev(SIZE_MAX), Constraint(SIZE_MAX), Orientation(true), Mark(false) {};

    Edge(size_t n, size_t s, Edge &nx, Edge &pr, Edge &tw) : Node(n), Self(s), Next(nx.Self), Twin(tw.Self), Prev(pr.Self), Constraint(SIZE_MAX), Orientation(true), Mark(false) {};

    Edge(size_t n, size_t s, size_t c, bool d) : Node(n), Self(s), Next(SIZE_MAX), Twin(SIZE_MAX), Prev(SIZE_MAX), Constraint(c), Orientation(d), Mark(false) {};

    size_t node() const { return Node; };

    size_t self() const { return Self; };

    size_t next() const { return Next; };

    size_t twin() const { return Twin; };

    size_t prev() const { return Prev; };

    bool orientation() const { return Orientation; };

    bool mark() const { return Mark; };

    bool operator==(Edge const &e) const {
        return (Node == e.Node &&
                Self == e.Self &&
                Next == e.Next &&
                Twin == e.Twin &&
                Prev == e.Prev &&
                Constraint == e.Constraint &&
                Orientation == e.Orientation);
    };

protected:
    size_t Node;            //Point at start of this edge
    size_t Self;            //This edge in this triangle
    size_t Next;            //Next edge in this triangle
    size_t Twin;            //Twin edge in adjacent triangle
    size_t Prev;            //Previous edge in this triangle

    size_t Constraint;      //DartConstraint index
    bool Orientation;       //don't care if unconstrained

    bool Mark;              //Auxiliary variable for mesh refinement
};

#endif //OERSTED_EDGE_H
