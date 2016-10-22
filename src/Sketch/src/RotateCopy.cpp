#include "Sketch.hpp"

RotateCopy::RotateCopy(std::vector<const Curve *> &input, Vertex *center, double angle, size_t copies) {
    // Creates rotated copies of the input curves about an vertex
    // #TODO: Need to rearrange code and reserve vector sizes in a way that makes more sense (much code copied from MirrorCopy constructor)
    // #TODO: Restructure to obviate the need for local_curves and local_verticies
    // #TODO: Three Groups: Leading Curves, Lagging Curves, Internal Curves
    // #TODO: Three Groups (?): Leading Verticies, Lagging Verticies, Internal Verticies

    // Assign properties
    Input = input;
    Center = center;
    Angle = angle;
    Copies = copies;
    Curves.reserve(copies * input.size());
    Verticies.reserve(3 * copies * input.size());
    Constraints.reserve(3 * copies * input.size());

    // Reserve Local Vectors
    std::vector<Curve *> local_curves;
    local_curves.reserve(Input.size());

    std::vector<Vertex *> local_verticies;
    local_verticies.reserve(Input.size() * 2);

    std::list<Vertex *> input_vlist;

    for (size_t i = 0; i != Copies; ++i) {
        local_curves.clear();
        local_verticies.clear();
        input_vlist.clear();

        // Clone input curves and create a list of unique input verticies
        for (auto c : Input) {
            local_curves.push_back(c->clone());
            Curves.push_back(local_curves.back());
            c->get_verticies(input_vlist);
        }
        input_vlist.sort();
        input_vlist.unique();

        // Rewrite
        /*
        // Get local curves
        for (auto c : Input) {
            local_curves.push_back(c->clone());
        }

        // Erase corotational curves
        for (size_t j = 0;j != local_curves.size();++j) {
            auto k = local_curves.begin();
            while (k != local_curves.end()) {
                if (local_curves[j].is_corotational(*k, center, angle)) {
                    k = local_curves.erase(k);
                }
                else {
                    ++k;
                }
            }
        }

        // Get local verticies
        for (auto c : local_curves) {
            c->get_verticies(input_vlist);
        }

        // Erase corotational verticies
        for (size_t j = 0;j != input_vlist.size();++j) {
            auto k = input_vlist.end();
            while (k != input_vlist.end()) {
                if (input_vlist[j].is_corotational(*k, center, angle)) {
                    k = input_vlist.erase(k);
                }
                else {
                    ++k;
                }
            }
        }
        */

        // Erase center if it exists
        auto j = std::find(input_vlist.begin(), input_vlist.end(), center);
        if (j != input_vlist.end()) {
            input_vlist.erase(j);
        }

        // Rotate input verticies
        double a = (i + 1.0) * Angle * M_PI / 180.0;
        double cosa = cos(a);
        double sina = sin(a);
        double x0 = Center->x();
        double y0 = Center->y();

        auto v = input_vlist.begin();
        while (v != input_vlist.end()) {
            double dx = (*v)->x() - x0;
            double dy = (*v)->y() - y0;

            double vx = cosa * dx - sina * dy + x0;
            double vy = sina * dx + cosa * dy + y0;

            local_verticies.push_back(new Vertex(vx, vy));

            Verticies.push_back(local_verticies.back());

            v++;
        }

        // Replace verticies in rotated curves
        std::vector<Vertex *> input_vvector{input_vlist.begin(), input_vlist.end()};
        for (auto c : local_curves) {
            c->replace_verticies(input_vvector, local_verticies);
        }

        // Constrain rotated verticies to be rotated about the origin
        for (size_t j = 0; j != local_verticies.size(); ++j) {
            Constraints.push_back(new Rotation(*input_vvector[j], *local_verticies[j], *Center, Angle * (i + 1)));
        }
    }

    // Shrink
    Curves.shrink_to_fit();
    Verticies.shrink_to_fit();
    Constraints.shrink_to_fit();
}