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

bool Edge::is_valid(Mesh const &mesh) const {
    bool value = true;

    value = value && (Self == mesh.Edges[Next]->Prev);
    value = value && (Self == mesh.Edges[Prev]->Next);
    value = value && (Self == mesh.Edges[Twin]->Twin);

    return value;
}

bool Edge::swap(Mesh const &mesh) {
    if (ConstraintCurve == nullptr) {
        Edge *e1 = mesh.Edges[Next];
        Edge *e2 = mesh.Edges[Prev];
        Edge *e3 = mesh.Edges[mesh.Edges[Twin]->Next];
        Edge *e4 = mesh.Edges[mesh.Edges[Twin]->Prev];
        Edge *twin = mesh.Edges[Twin];

        Node = e2->Node;
        Next = e4->Self;
        Prev = e1->Self;
        Mark = false;

        twin->Node = e4->Node;
        twin->Next = e2->Self;
        twin->Prev = e3->Self;
        twin->Mark = false;

        e1->Next = Self;
        e1->Prev = e4->Self;
        e1->Mark = false;

        e2->Next = e3->Self;
        e2->Prev = Twin;
        e2->Mark = false;

        e3->Next = Twin;
        e3->Prev = e2->Self;
        e3->Mark = false;

        e4->Next = e1->Self;
        e4->Prev = Self;
        e4->Mark = false;

        return true;
    } else {
        Mark = false;
        return false;
    }
}

void Edge::recursive_mark(Mesh const &mesh) {
    if (Mark && mesh.Edges[Next]->Mark && mesh.Edges[Prev]->Mark) {
        mesh.Edges[Next]->Mark = false;
        mesh.Edges[Prev]->Mark = false;
    }
}