#include "Constellation.h"
#include "Sketch.h"
#include "Star.h"
#include "Branch.h"
#include "Curve.h"
#include "Contour.h"

Constellation::Constellation(Sketch const *s) {
    for (size_t i = 0; i != s->size_verticies(); ++i) {
        Stars.push_back(Star{s->vertex(i), s});
    }
    pop();
}

bool Constellation::twin(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out) {
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
    std::pair<double,double> sup{0.0,0.0};
    double ang{0.0};

    for (auto s = Stars.begin(); s != Stars.end(); ++s) {
        for (auto b = s->begin(); b != s->end(); ++b) {
            std::pair<double,double> sup_ij = b->Path->supremum();
            if ((sup_ij > sup) || (sup == sup_ij && b->Angle < ang)) {
                sup = sup_ij;
                ang = b->Angle;

                s_out = s;
                b_out = b;
            }
        }
    }
}

void Constellation::pop(std::shared_ptr<Curve> c) {
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

std::shared_ptr<Contour> Constellation::boundary() {
    std::vector<std::shared_ptr<Curve>> curves;
    std::vector<bool> orientation;

    // Base of induction
    auto s{Stars.begin()};
    auto b{s->begin()};

    supremum(s,b);

    double angle{0.0};
    b = s->prev(b);

    curves.push_back(b->Path);
    orientation.push_back(b->Orientation);

    while (true) {
        bool has_twin = twin(s, b);

        if (!has_twin) {
            return std::make_shared<Contour>();
        } else {
            b = s->prev(b);
            angle += 2.0 * M_PI - b->Angle;
        }

        if (b->Path == curves.back()) { //excludes self-closed contours
            return std::make_shared<Contour>();
        } else if (b->Path == curves[0]) {
            double tol = (FLT_EPSILON * M_PI * (curves.size() - 1));
            double expected = (curves.size() - 2) * M_PI;
            double diff = abs(angle - expected);

            if (diff <= tol) {
                return std::make_shared<Contour>(curves,orientation);
            } else {
                return std::make_shared<Contour>();
            }
        } else {
            curves.push_back(b->Path);
            orientation.push_back(b->Orientation);
        }
    }
}

std::vector<std::shared_ptr<Contour>> Constellation::contours() {
    std::vector<std::shared_ptr<Contour>> contours;
    std::vector<std::shared_ptr<Curve>> contour_curves;
    std::vector<bool> orientation;

    while (size() > 0) {
        if (find_closed_contour(contour_curves, orientation)) {
            contours.push_back(std::make_shared<Contour>(contour_curves, orientation));
        } else {
            return std::vector<std::shared_ptr<Contour>>(); // TODO: Ugly multiple return points
        }
    }

    return contours;
}

bool Constellation::find_closed_contour(std::vector<std::shared_ptr<Curve>> &curves, std::vector<bool> &orientation) {
    curves.resize(0);
    orientation.resize(0);

    // Base of induction
    auto s{Stars.begin()};
    auto b{s->begin()};

    supremum(s,b);

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