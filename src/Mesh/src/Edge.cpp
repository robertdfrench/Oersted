#include "Mesh.hpp"

Edge::Edge(Curve *c, bool direction, size_t v) {
    ConstraintCurve = c;
    Node = v;
    Twin = this;
    Next = nullptr;
    Prev = nullptr;
    Orientation = direction;
}

bool Edge::is_valid() const {
    bool value = true;

    value = value && (this == Next->Prev);
    value = value && (this == Prev->Next);
    value = value && (this == Twin->Twin);

    return value;
}

bool Edge::swap() {
    if (ConstraintCurve == nullptr) {
        Edge * e1 = Next;
        Edge * e2 = Prev;
        Edge * e3 = Twin->Next;
        Edge * e4 = Twin->Prev;

        Node = e2->Node;
        Next = e4;
        Prev = e1;
        Mark = false;

        Twin->Node = e4->Node;
        Twin->Next = e2;
        Twin->Prev = e3;
        Twin->Mark = false;

        e1->Next = this;
        e1->Prev = e4;
        e1->Mark = false;

        e2->Next = e3;
        e2->Prev = Twin;
        e2->Mark = false;

        e3->Next = Twin;
        e3->Prev = e2;
        e3->Mark = false;

        e4->Next = e1;
        e4->Prev = this;
        e4->Mark = false;

        return true;
    } else {
        Mark = false;
        return false;
    }
}

void Edge::recursive_mark() {
    if (Mark && Next->Mark && Prev->Mark) {
        Next->Mark = false;
        Prev->Mark = false;
    }
}