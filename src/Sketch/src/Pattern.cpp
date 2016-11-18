#include "Sketch.hpp"

void Pattern::register_elements(Sketch *s) {
    for (auto v : Verticies) {
        s->add_element(v);
    }

    for (auto c : Curves) {
        s->add_element(*c);
    }

    for (auto c : Constraints) {
        s->add_element(*c);
    }
}