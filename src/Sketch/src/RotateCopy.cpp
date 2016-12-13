#include "RotateCopy.h"
#include "Curve.h"
#include "Rotation.h"

RotateCopy::RotateCopy(std::vector<std::shared_ptr<Curve const>> input, std::shared_ptr<Vertex const> center, double angle, size_t copies, bool remove_internal) {
    // Creates rotated copies of the input curves about an vertex
    // TODO: Check for complete elimination of leading/lagging curves

    // Assign properties
    Input = input;
    Center = center;
    Angle = angle;
    Copies = copies;
    Curves.reserve(copies * input.size());
    Verticies.reserve(3 * copies * input.size());
    Constraints.reserve(3 * copies * input.size());
    RemoveInternalBoundaries = remove_internal;

    std::vector<std::shared_ptr<Curve const>> leading_curves;
    std::vector<std::shared_ptr<Curve const>> lagging_curves;
    std::vector<bool> is_lead_or_lag(input.size(), false);

    // Find leading/lagging curve and vertex pairs
    std::list<std::shared_ptr<Vertex const>> lead_vlist;
    std::list<std::shared_ptr<Vertex const>> lag_vlist;
    for (size_t i = 0; i != input.size(); ++i) {
        for (size_t j = 0; j != input.size(); ++j) {
            if (!(i == j || is_lead_or_lag[i] || is_lead_or_lag[j])) {
                MatchOrientation dir = input[i]->is_identical(input[j], center, angle);
                if(dir != MatchOrientation::None){
                    leading_curves.push_back(input[i]);
                    is_lead_or_lag[i] = true;
                    input[i]->get_verticies(lead_vlist);

                    // dir indicates order which lagging curve verticies match the leading curve verticies
                    lagging_curves.push_back(input[j]);
                    is_lead_or_lag[j] = true;
                    input[j]->get_verticies(lag_vlist, dir);
                }
            }
        }
    }

    // Sort leading/lagging verticies
    std::vector<std::shared_ptr<Vertex const>> leading_verticies;
    std::vector<std::shared_ptr<Vertex const>> lagging_verticies;
    {
        lead_vlist.remove(center);
        lag_vlist.remove(center);

        std::list<std::pair<std::shared_ptr<Vertex const>, std::shared_ptr<Vertex const>>> vpairs;
        auto j = lag_vlist.begin();
        for (auto i : lead_vlist) {
            vpairs.push_back(std::pair<std::shared_ptr<Vertex const>, std::shared_ptr<Vertex const>>(i, *j));
            ++j;
        }
        vpairs.sort();
        vpairs.unique();

        for (auto &i : vpairs) {
            leading_verticies.push_back(i.first);
            lagging_verticies.push_back(i.second);
        }
    }

    // If not leading/lagging, curve is internal
    std::vector<std::shared_ptr<Curve const>> internal_curves;
    for (size_t i = 0; i != input.size(); ++i) {
        if (!is_lead_or_lag[i]) {
            internal_curves.push_back(input[i]);
        }
    }

    // Get internal verticies
    std::vector<std::shared_ptr<Vertex const>> internal_verticies;
    {
        std::list<std::shared_ptr<Vertex const>> local_vlist;
        for (auto c : internal_curves) {
            c->get_verticies(local_vlist);
        }
        local_vlist.sort(); // remove duplicates
        local_vlist.unique();
        local_vlist.remove(center);

        // Trim leading verticies
        for (auto v : leading_verticies) {
            local_vlist.remove(v);
        }

        // Trim lagging verticies
        for (auto v : lagging_verticies) {
            local_vlist.remove(v);
        }

        // Copy to vector
        internal_verticies.assign(local_vlist.begin(), local_vlist.end());
    }

    // Make rotated copies
    std::vector<std::shared_ptr<Vertex const>> rotated_lagging(leading_verticies.begin(), leading_verticies.end());
    for (size_t i = 0; i != Copies; ++i) {
        bool last_iteration = (i == Copies - 1);
        // Create curve clones
        std::vector<std::shared_ptr<Curve>> local_curves;
        for (auto c : internal_curves) {
            local_curves.push_back(c->clone());
            Curves.push_back(local_curves.back());
        }

        for (auto c : leading_curves) {
            local_curves.push_back(c->clone());
            Curves.push_back(local_curves.back());

            if (!last_iteration) {
                std::const_pointer_cast<Curve>(Curves.back())->for_construction(RemoveInternalBoundaries);
            } else {
                std::const_pointer_cast<Curve>(c)->for_construction(RemoveInternalBoundaries);
            }
        }

        // Create new verticies and constrain them
        double a = (i + 1) * Angle * M_PI / 180.0;
        double cosa = cos(a);
        double sina = sin(a);
        double x0 = Center->x();
        double y0 = Center->y();

        std::vector<std::shared_ptr<Vertex const>> rotated_leading;
        for (auto v : leading_verticies) {
            double dx = v->x() - x0;
            double dy = v->y() - y0;

            double vx = cosa * dx - sina * dy + x0;
            double vy = sina * dx + cosa * dy + y0;

            Verticies.push_back(std::make_shared<Vertex>(vx, vy));

            rotated_leading.push_back(Verticies.back());

            Constraints.push_back(std::make_shared<Rotation>(v, rotated_leading.back(), Center, Angle * (i + 1)));
        }

        std::vector<std::shared_ptr<Vertex const>> rotated_internal;
        for (auto v : internal_verticies) {
            double dx = v->x() - x0;
            double dy = v->y() - y0;

            double vx = cosa * dx - sina * dy + x0;
            double vy = sina * dx + cosa * dy + y0;

            Verticies.push_back(std::make_shared<Vertex>(vx, vy));

            rotated_internal.push_back(Verticies.back());

            Constraints.push_back(std::make_shared<Rotation>(v, rotated_internal.back(), Center, Angle * (i + 1)));
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