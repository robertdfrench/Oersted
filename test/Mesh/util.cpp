#include "test_Mesh.hpp"

bool edges_are_optimal(Mesh &m) {
    for (size_t i = 0;i < m.size_edges();++i) {
        EXPECT_TRUE(m.edge(i)->is_locally_optimal());
    }

    return true;
}

bool edges_are_valid(Mesh &m) {
    for (size_t i = 0;i < m.size_edges();++i) {
        const Edge *e = m.edge(i);

        EXPECT_TRUE(e == e->next()->prev());
        EXPECT_TRUE(e == e->prev()->next());
        EXPECT_TRUE(e == e->twin()->twin());

        //e->ConstraintCurve == e->Twin->ConstraintCurve and e->Orientation != e->Orientation
        if (!(e->twin() == e)) {
            EXPECT_TRUE(e->node() == e->twin()->next()->node());
            EXPECT_TRUE(e->constraint_curve() == e->twin()->constraint_curve());
            if (e->constraint_curve() != nullptr) {
                EXPECT_TRUE(e->orientation() != e->twin()->orientation());
            }

            EXPECT_FALSE(e->node() == e->twin()->node());
        }

        if (e->constraint_curve() != nullptr) {
            if (e->orientation()) {
                EXPECT_TRUE(*e->base() == *e->constraint_curve()->start());
                EXPECT_TRUE(*e->tip() == *e->constraint_curve()->end());
            }
            else {
                EXPECT_TRUE(*e->base() == *e->constraint_curve()->end());
                EXPECT_TRUE(*e->tip() == *e->constraint_curve()->start());
            }
        }
    }

    return true;
}

void forced_refinement(Mesh &m, std::string file_name, size_t num_refines) {
    double minq = m.MinimumElementQuality;
    double mine = m.MinimumElementSize;
    double maxe = m.MaximumElementSize;

    m.MinimumElementQuality = DBL_MAX;
    m.MinimumElementSize = 0.0;
    m.MaximumElementSize = DBL_MAX;

    for (size_t i = 0;i < num_refines;++i) {
        try {
            m.refine_once();
        }
        catch (const std::exception except) {
            throw;
        }
        m.save_as(SAVE_DIR, file_name + "_" + std::to_string(i));

        edges_are_valid(m);
        edges_are_optimal(m);
    }

    m.MinimumElementQuality = minq;
    m.MaximumElementSize = mine;
    m.MaximumElementSize = maxe;
}

std::vector<size_t> map_verticies_to_points(std::vector<Vertex> verts, Mesh m) {
    std::vector<size_t> index;
    index.reserve(verts.size());

    for (size_t i = 0;i != verts.size();++i) {
        size_t j = 0;
        for (size_t j = 0;j != m.size_points();++j) {
            if (verts[i].x() == m.point(j)->X && verts[i].y() == m.point(j)->Y) {
                index.push_back(j);
                break;
            }
        }
    }

    return std::move(index);
}