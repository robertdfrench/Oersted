#include "Sketch.hpp"

Constellation::Constellation(const Sketch *s) {
    for (size_t i = 0; i != s->size_verticies(); ++i) {
        Stars.push_back(Star{s->vertex(i), s});
    }
    pop();
}

bool Constellation::twin(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out) {
    const Curve *c = b_out->Path; // #TODO: Local variable 'c' is only assigned but never accessed

    for (auto s = Stars.begin(); s != Stars.end(); ++s) {
        if (s != s_out) {
            for (auto b = s->begin(); b != s->end(); ++b) {
                if (b->Path == b_out->Path) {
                    b_out = b;
                    s_out = s;

                    return true;
                }
            }
        }
    }

    return false;
}

void Constellation::supremum(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out) {
    double sup = 0.0;
    double ang;

    for (auto s = Stars.begin(); s != Stars.end(); ++s) {
        for (auto b = s->begin(); b != s->end(); ++b) {
            double sup_ij = b->Path->supremum();
            if (sup < sup_ij || (sup == sup_ij && b->Angle < ang)) {
                sup = sup_ij;
                ang = b->Angle;

                s_out = s;
                b_out = b;
            }
        }
    }
}

void Constellation::pop(const Curve *c) {
    bool iterate = true;

    while (iterate) {
        iterate = false;

        // remove curve c
        for (auto s = Stars.begin(); s != Stars.end(); ++s) {
            s->pop(c);
        }

        auto s = Stars.begin();
        while (s != Stars.end()) {
            if (s->size() == 0) {
                s = Stars.erase(s);
            } else {
                ++s;
            }
        }

        // pop curves from stars which only have one branch
        for (auto s = Stars.begin(); s != Stars.end(); ++s) {
            if (s->size() == 1) {
                c = s->begin()->Path;
                iterate = true;
                break;
            }
        }
    }
}

bool Constellation::boundary(Contour *c) {
    std::vector<const Curve *> curves;
    std::vector<bool> orientation;

    // Base of induction
    auto s{Stars.begin()};
    auto b{s->begin()};

    supremum(s, b);

    double angle{0.0};
    b = s->prev(b);

    curves.push_back(b->Path);
    orientation.push_back(b->Orientation);

    while (true) {
        bool has_twin = twin(s, b);

        if (!has_twin) {
            return false;
        } else {
            b = s->prev(b);
            angle += 2.0 * M_PI - b->Angle;
        }

        if (b->Path == curves.back()) { //excludes self-closed contours
            return false;
        } else if (b->Path == curves[0]) {
            double tol = (FLT_EPSILON * M_PI * (curves.size() - 1));
            double expected = (curves.size() - 2) * M_PI;
            double diff = abs(angle - expected);

            if (diff <= tol) {
                c->initialize(curves, orientation);
                return true;
            } else {
                return false;
            }
        } else {
            curves.push_back(b->Path);
            orientation.push_back(b->Orientation);
        }
    }
}

bool Constellation::contours(std::vector<Contour *> &contours) {
    std::vector<const Curve *> contour_curves;
    std::vector<bool> orientation;

    while (size() > 0) {
        bool success = find_closed_contour(contour_curves, orientation);
        if (success) {
            contours.push_back(new Contour(contour_curves, orientation));
        } else {
            return false;
        }
    }

    return true;
}

bool Constellation::find_closed_contour(std::vector<const Curve *> &curves, std::vector<bool> &orientation) {
    curves.resize(0);
    orientation.resize(0);

    // Base of induction
    auto s{Stars.begin()};
    auto b{s->begin()};

    supremum(s, b);

    double angle{0.0};
    b = s->next(b);

    curves.push_back(b->Path);
    orientation.push_back(b->Orientation);

    while (true) {
        bool has_twin = twin(s, b);

        if (!has_twin) {
            return false;
        } else {
            angle += b->Angle;
            b = s->next(b);
        }

        if (b->Path == curves.back()) { //excludes self-closed contours
            return false;
        } else if (b->Path == curves[0]) {
            double tol = (FLT_EPSILON * M_PI * (curves.size() - 1));
            double expected = (curves.size() - 2) * M_PI;
            double diff = abs(angle - expected);

            if (diff <= tol) {
                pop(curves.back());
                return true;
            } else {
                return false;
            }
        } else {
            curves.push_back(b->Path);
            orientation.push_back(b->Orientation);
        }
    }
}