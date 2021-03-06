#include "MirrorCopy.h"
#include "Vertex.h"
#include "LineSegment.h"
#include "Symmetry.h"

MirrorCopy::MirrorCopy(std::vector<std::shared_ptr<Curve const>> &input, std::shared_ptr<LineSegment const> l, bool remove_internal) {
    // Creates mirror copies of the input curves about a line

    // Assign Properties
    Input = input;
    SymmetryLine = l;
    RemoveInternalBoundaries = remove_internal;

    // Clone input curves and create a list of unique input verticies
    Curves.reserve(Input.size());

    std::list<std::shared_ptr<Vertex const>> input_vlist;
    for (auto c : Input) {
        if (l->is_coincident(c)) {
            std::const_pointer_cast<Curve>(c)->for_construction(RemoveInternalBoundaries);
        } else {
            Curves.push_back(c->clone());
            c->get_verticies(input_vlist);
        }
    }
    input_vlist.sort();
    input_vlist.unique();

    // Mirror input verticies
    double x0 = l->start()->x();
    double y0 = l->start()->y();
    double dx = l->end()->x() - x0;
    double dy = l->end()->y() - y0;
    double d = sqrt(dx * dx + dy * dy);
    dx /= d;
    dy /= d;

    Verticies.reserve(input_vlist.size());
    auto v = input_vlist.begin();
    while (v != input_vlist.end()) {
        double x = (*v)->x();
        double y = (*v)->y();
        double px = x - x0;
        double py = y - y0;
        double pd = dx * px + dy * py;

        px -= pd * dx;
        py -= pd * dy;

        if (sqrt(px * px + py * py) > d * FLT_EPSILON) {
            px = x - 2.0 * px;
            py = y - 2.0 * py;

            //Verticies.push_back(new Vertex(px, py));
            Verticies.push_back(std::make_shared<Vertex>(px, py));


            ++v;
        } else {
            v = input_vlist.erase(v);
        }
    }

    // Replace verticies in mirror curves
    std::vector<std::shared_ptr<Vertex const>> input_vvector{input_vlist.begin(), input_vlist.end()};
    for (auto c : Curves) {
        c->replace_verticies(input_vvector, std::vector<std::shared_ptr<Vertex const>>{Verticies.begin(), Verticies.end()});
        c->reverse();
    }

    // Constrain mirrored verticies to be symmetric about the SymmetryLine
    Constraints.reserve(Verticies.size());
    for (size_t i = 0; i != Verticies.size(); ++i) {
        Constraints.push_back(std::make_shared<Symmetry>(input_vvector[i], Verticies[i], SymmetryLine));
    }
}