#include "test_Sketch.hpp"

void test_sketch_size(Sketch &s, size_t nverts, size_t ncurves, size_t nconstraints, size_t ncontours) {
    EXPECT_EQ(nverts, s.size_verticies());
    EXPECT_EQ(ncurves, s.size_curves());
    EXPECT_EQ(nconstraints, s.size_constraints());
    EXPECT_EQ(ncontours, s.size_contours());
}

bool has_rotational_image(Sketch &s, std::shared_ptr<Vertex const> const v, std::shared_ptr<Vertex const> const center, double angle) {
    double2 p = v->rotate(center, angle);

    bool rotation_found = false;
    for (size_t j = 0; j != s.size_verticies(); ++j) {
        std::shared_ptr<Vertex const> vj = s.vertex(j);
        rotation_found = (abs(vj->x() - p.X) < TOL && abs(vj->y() - p.Y) < TOL);

        if (rotation_found) {
            break;
        }
    }

    return rotation_found;
}

bool has_rotational_image(Sketch &s, std::shared_ptr<Vertex const> const v0, std::shared_ptr<Vertex const> const v1, std::shared_ptr<Vertex const> const center, double angle) {
    double2 p0 = v0->rotate(center, angle);

    double2 p1 = v1->rotate(center, angle);

    bool rotation0_found = false;
    bool rotation1_found = false;

    for (size_t j = 0; j != s.size_curves(); ++j) {
        auto c = s.curve(j);
        auto vs = c->start();
        auto ve = c->end();

        rotation0_found = (abs(vs->x() - p0.X) < TOL && abs(vs->y() - p0.Y) < TOL);

        if (rotation0_found) {
            rotation1_found = (abs(ve->x() - p1.X) < TOL && abs(ve->y() - p1.Y) < TOL);
        } else {
            rotation0_found = (abs(ve->x() - p0.X) < TOL && abs(ve->y() - p0.Y) < TOL);
            if (rotation0_found) {
                rotation1_found = (abs(vs->x() - p1.X) < TOL && abs(vs->y() - p1.Y) < TOL);
            }
        }

        if (rotation0_found && rotation1_found) {
            break;
        }
    }

    return rotation0_found && rotation1_found;
}

void test_rotated_verticies(Sketch &s, std::vector<size_t> index, std::shared_ptr<Vertex const> const center, double angle, size_t copies) {
    for (size_t j = 1; j != (copies + 1); ++j) {
        for (size_t i : index) {
            EXPECT_TRUE(has_rotational_image(s, s.vertex(i), center, angle * j));
        }
    }
}

void test_rotated_curves(Sketch &s, std::vector<size_t> index, std::shared_ptr<Vertex const> const center, double angle, size_t copies) {
    for (size_t j = 1; j != (copies + 1); ++j) {
        for (size_t i : index) {
            EXPECT_TRUE(has_rotational_image(s, s.curve(i)->start(), s.curve(i)->end(), center, angle * j));
        }
    }
}


TEST(Sketch, constructor) {
    Sketch s;
}

TEST(Sketch, new_element) {
    { // Case 0: Vertex
        Sketch s;

        auto v0 = s.new_element<Vertex>();
        EXPECT_TRUE(s.size() == 1);
        EXPECT_TRUE(v0->get_equation_index() == 0);
        EXPECT_TRUE(v0->x_index() == 0);
        EXPECT_TRUE(v0->y_index() == 1);

        auto v1 = s.new_element<Vertex>();

        EXPECT_TRUE(s.size() == 2);
        EXPECT_TRUE(v1->get_equation_index() == 0);
        EXPECT_TRUE(v1->x_index() == 2);
        EXPECT_TRUE(v1->y_index() == 3);
    }

    { // Case 1: LineSegment
        Sketch s;

        auto line = s.new_element<LineSegment>();
        EXPECT_TRUE(s.size() == 1);
    }

    { // Case 2: CircularArc
        Sketch s;

        auto line = s.new_element<LineSegment>();
        EXPECT_TRUE(s.size() == 1);
    }

    { // Case 3: Length
        Sketch s;

        auto v0 = s.new_element<Vertex>(3.14159, 2.7183);
        auto v1 = s.new_element<Vertex>(6.14159, 6.7183);

        auto line = s.new_element<LineSegment>(v0, v1);

        auto length = s.new_element<Length>(line, 1.0);

        EXPECT_TRUE(s.size() == 4);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Sketch__new_element_Length");
    }
}