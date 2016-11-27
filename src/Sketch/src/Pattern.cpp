#include "Pattern.h"
#include "Sketch.h"
#include "Vertex.h"
#include "Curve.h"
#include "Constraint.h"

void Pattern::register_elements(Sketch *s) const {
    for (auto v : Verticies) {
        s->add_element(v);
    }

    for (auto c : Curves) {
        s->add_element(c);
    }

    for (auto c : Constraints) {
        s->add_element(c);
    }
}