#include <numeric>

#include "Branch.h"
#include "Curve.h"
#include "Star.h"
#include "Sketch.h"
#include "Vertex.h"

Star::Star(std::shared_ptr<Vertex> v, Sketch const *s) {
    StarVertex = v;

    // Extract Curves
    std::vector<std::shared_ptr<Curve>> curves;
    std::vector<bool> orientation;
    for (size_t i = 0; i < s->size_curves(); ++i) {
        std::shared_ptr<Curve> c = s->curve(i);
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

    double tol = M_PI * FLT_EPSILON;
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
        std::shared_ptr<Vertex> v;
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

        Branches.push_back(Branch{curves[j], angle[j] - angle[k], orientation[j]});

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

std::shared_ptr<Curve> Star::next(std::shared_ptr<Curve> c) const {
    for (auto b = begin(); b != end(); ++b) {
        if (b->Path == c) {
            return next(b)->Path;
        }
    }

    return nullptr;
}

std::shared_ptr<Curve> Star::prev(std::shared_ptr<Curve> c) const {
    for (auto b = begin(); b != end(); ++b) {
        if (b->Path == c) {
            return prev(b)->Path;
        }
    }

    return nullptr;
}

void Star::pop(std::shared_ptr<Curve> c) {
    for (auto b = begin(); b != end(); ++b) {
        if (b->Path == c) {
            prev(b)->Angle += b->Angle;

            Branches.erase(b);

            break;
        }
    }
}