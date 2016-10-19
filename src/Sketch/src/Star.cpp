#include "Sketch.hpp"

Star::Star(const Vertex *v, const Sketch *s) {
    StarVertex = v;

    // Extract Curves
    std::vector<const Curve *> curves;
    std::vector<bool> orientation;
    for (size_t i = 0; i < s->size_curves(); ++i) {
        const Curve *c = s->curve(i);
        if (!c->ForConstruction) {
            if (StarVertex == c->start()) {
                curves.push_back(c);
                orientation.push_back(true);
            } else if (StarVertex == c->end()) {
                curves.push_back(c);
                orientation.push_back(false);
            }
        }
    }

    // Calculate angle tangent to star point
    std::vector<double> angle;
    std::vector<double> d_angle;
    for (size_t i = 0; i < curves.size(); ++i) {
        if (orientation[i]) {
            angle.push_back(curves[i]->a(0.0, true));
            d_angle.push_back(curves[i]->da(0.0, true));
        } else {
            angle.push_back(curves[i]->a(1.0, false));
            d_angle.push_back(curves[i]->da(1.0, false));
        }
    }

    // Sort by tangent angle
    std::vector<size_t> index;
    index.resize(angle.size());
    std::iota(index.begin(), index.end(), 0);

    double tol = M_PI * sqrt(DBL_EPSILON);
    std::sort(index.begin(), index.end(),
              [&](size_t i, size_t j) {
                  if (abs(angle[i] - angle[j]) < tol ||
                      (abs(abs(angle[i]) - M_PI) < tol && abs(abs(angle[j]) - M_PI) < tol)) {
                      return (d_angle[i] > d_angle[j]);
                  } else {
                      return (angle[i] > angle[j]);
                  }
              }
    );

    // Store end point angle difference
    double x = StarVertex->x();
    double y = StarVertex->y();
    for (size_t i = 0; i != curves.size(); ++i) {
        const Vertex *v;
        if (orientation[i]) {
            v = curves[i]->end();
        } else {
            v = curves[i]->start();
        }

        double dx = v->x() - x;
        double dy = v->y() - y;

        angle[i] = atan2(dy, dx);
    }

    for (size_t i = 0; i != index.size(); ++i) {
        size_t j = index[i];
        size_t k = index[(i + 1) % index.size()];

        Branches.push_back(Branch{curves[j], orientation[j], angle[j] - angle[k]});

        // Handle branch cut
        if (Branches.back().Angle < 0.0) {
            Branches.back().Angle += 2.0 * M_PI;
        }
    }

    // Handle self-closed Star with two branches having shared endpoints
    if (index.size() == 2 && Branches.back().Angle == 0.0) {
        Branches.back().Angle += 2.0 * M_PI;
    }

}

const Curve *Star::next(const Curve *c) const {
    for (auto b = begin(); b != end(); ++b) {
        if (b->Path == c) {
            return next(b)->Path;
        }
    }

    return nullptr;
}

const Curve *Star::prev(const Curve *c) const {
    for (auto b = begin(); b != end(); ++b) {
        if (b->Path == c) {
            return prev(b)->Path;
        }
    }

    return nullptr;
}

void Star::pop(const Curve *c) {
    for (auto b = begin(); b != end(); ++b) {
        if (b->Path == c) {
            prev(b)->Angle += b->Angle;

            Branches.erase(b);

            break;
        }
    }
}