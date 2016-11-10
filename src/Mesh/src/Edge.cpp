#include "Mesh.hpp"

Edge::Edge(Curve *c, bool d, size_t n, size_t s) {
    ConstraintCurve = c;
    Node = n;
    Self = s;
    Twin = SIZE_MAX;
    Next = SIZE_MAX;
    Prev = SIZE_MAX;
    Orientation = d;
}

size_t Edge::tip(Mesh const &mesh) const {
    return (Next == Self ? mesh.Edges[Twin]->Node : mesh.Edges[Next]->Node);
};