#include "Sketch.hpp"

RotateCopy::RotateCopy(std::vector<const Curve *> &input, Vertex *center, double angle, size_t copies, bool remove_internal) {
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
    RemoveInternalBoundaries = remove_internal;

    std::vector<const Curve *> leading_curves;
    std::vector<const Curve *> lagging_curves;
    std::vector<bool> is_lead_or_lag(input.size(),false);

    // Find leading/lagging curve pairs
    for (size_t i = 0; i != input.size(); ++i) {
        for (size_t j = 0; j != input.size(); ++j) {
            if (!(i == j || is_lead_or_lag[i] || is_lead_or_lag[j])) { //
                if (input[i]->is_identical(input[j], center, angle)) {
                    leading_curves.push_back(input[i]);
                    is_lead_or_lag[i] = true;

                    lagging_curves.push_back(input[j]);
                    is_lead_or_lag[j] = true;
                    break;
                }
            }
        }
    }

    // If not leading/lagging, curve is internal
    std::vector<const Curve *> internal_curves;
    for (size_t i =0; i != input.size(); ++i) {
        if(!is_lead_or_lag[i]) {
            internal_curves.push_back(input[i]);
        }
    }

    // TODO: Check for complete elimination of leading/lagging curves

    // Get leading/lagging verticies
    std::vector<Vertex *> leading_verticies;
    {
        std::list<Vertex *> local_verts;
        for (auto i : leading_curves) {
            i->get_verticies(local_verts);
        }
        local_verts.sort(); // remove duplicates
        local_verts.unique();
        local_verts.remove(center);

        leading_verticies.assign(local_verts.begin(),local_verts.end()); // assign to vector
    }

    std::vector<Vertex *> lagging_verticies;
    {
        std::list<Vertex *> local_verts;
        for (auto i : lagging_curves) {
            i->get_verticies(local_verts);
        }
        local_verts.sort(); // remove duplicates
        local_verts.unique();
        local_verts.remove(center);

        lagging_verticies.assign(local_verts.begin(),local_verts.end()); // assign to vector
    }

    // Get internal verticies
    std::vector<Vertex *> internal_verticies;
    {
        std::list<Vertex *> local_verts;
        for (auto c : internal_curves) {
            c->get_verticies(local_verts);
        }
        local_verts.sort(); // remove duplicates
        local_verts.unique();
        local_verts.remove(center);

        // Trim leading verticies
        for(auto v : leading_verticies) {
            local_verts.remove(v);
        }

        // Trim lagging verticies
        for(auto v : lagging_verticies) {
            local_verts.remove(v);
        }

        // Copy to vector
        internal_verticies.assign(local_verts.begin(),local_verts.end());
    }

    // Make rotated copies
    std::vector<Vertex *> rotated_lagging(leading_verticies.begin(), leading_verticies.end());
    for (size_t i = 0; i != Copies; ++i) {
        bool last_iteration = (i == Copies - 1);
        // Create curve clones
        std::vector<Curve *> local_curves;
        for (auto c : internal_curves) {
            local_curves.push_back(c->clone());
            Curves.push_back(local_curves.back());
        }

        for (auto c : leading_curves) {
            local_curves.push_back(c->clone());
            Curves.push_back(local_curves.back());

            if (RemoveInternalBoundaries) {
                if (!last_iteration) {
                    Curves.back()->ForConstruction = true;
                } else {
                    const_cast<Curve *>(c)->ForConstruction = true; // TODO: const_cast is ugly
                }
            }
        }

        // Create new verticies and constrain them
        double a = (i + 1) * Angle * M_PI / 180.0;
        double cosa = cos(a);
        double sina = sin(a);
        double x0 = Center->x();
        double y0 = Center->y();

        std::vector<Vertex *> rotated_leading;
        for(auto v : leading_verticies) {
            double dx = v->x() - x0;
            double dy = v->y() - y0;

            double vx = cosa * dx - sina * dy + x0;
            double vy = sina * dx + cosa * dy + y0;

            rotated_leading.push_back(new Vertex(vx, vy));

            Verticies.push_back(rotated_leading.back());

            Constraints.push_back(new Rotation(*v, *rotated_leading.back(), *Center, Angle * (i + 1)));
        }

        std::vector<Vertex *> rotated_internal;
        for(auto v : internal_verticies) {
            double dx = v->x() - x0;
            double dy = v->y() - y0;

            double vx = cosa * dx - sina * dy + x0;
            double vy = sina * dx + cosa * dy + y0;

            rotated_internal.push_back(new Vertex(vx, vy));

            Verticies.push_back(rotated_internal.back());

            Constraints.push_back(new Rotation(*v, *rotated_internal.back(), *Center, Angle * (i + 1)));
        }

        // Replace verticies in clones curves
        for (auto c : local_curves) {
            c->replace_verticies(leading_verticies, rotated_leading);
            c->replace_verticies(internal_verticies, rotated_internal);
            c->replace_verticies(lagging_verticies, rotated_lagging);
        }

        // Leading verticies becomes lagging verticies
        rotated_lagging.assign(rotated_leading.begin(), rotated_leading.end());
    }

    Curves.shrink_to_fit();
    Verticies.shrink_to_fit();
    Constraints.shrink_to_fit();
}