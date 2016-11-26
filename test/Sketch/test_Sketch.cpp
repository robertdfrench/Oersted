#include "test_Sketch.hpp"

void test_sketch_size(Sketch &s, size_t nverts, size_t ncurves, size_t nconstraints, size_t ncontours) {
    EXPECT_EQ(nverts, s.size_verticies());
    EXPECT_EQ(ncurves, s.size_curves());
    EXPECT_EQ(nconstraints, s.size_constraints());
    EXPECT_EQ(ncontours, s.size_contours());
}

bool has_rotational_image(Sketch &s, std::shared_ptr<Vertex> v, std::shared_ptr<Vertex> center, double angle) {
    double x, y;
    std::tie(x, y) = v->rotate(center, angle);

    bool rotation_found = false;
    for (size_t j = 0; j != s.size_verticies(); ++j) {
        std::shared_ptr<Vertex> vj = s.vertex(j);
        rotation_found = (abs(vj->x() - x) < TOL && abs(vj->y() - y) < TOL);

        if (rotation_found) {
            break;
        }
    }

    return rotation_found;
}

bool has_rotational_image(Sketch &s, std::shared_ptr<Vertex> v0, std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> center, double angle) {
    double x0, y0;
    std::tie(x0, y0) = v0->rotate(center, angle);

    double x1, y1;
    std::tie(x1, y1) = v1->rotate(center, angle);

    bool rotation0_found = false;
    bool rotation1_found = false;

    for (size_t j = 0; j != s.size_curves(); ++j) {
        std::shared_ptr<Curve> c = s.curve(j);
        std::shared_ptr<Vertex> vs = c->start();
        std::shared_ptr<Vertex> ve = c->end();

        rotation0_found = (abs(vs->x() - x0) < TOL && abs(vs->y() - y0) < TOL);

        if (rotation0_found) {
            rotation1_found = (abs(ve->x() - x1) < TOL && abs(ve->y() - y1) < TOL);
        } else {
            rotation0_found = (abs(ve->x() - x0) < TOL && abs(ve->y() - y0) < TOL);
            if (rotation0_found) {
                rotation1_found = (abs(vs->x() - x1) < TOL && abs(vs->y() - y1) < TOL);
            }
        }

        if (rotation0_found && rotation1_found) {
            break;
        }
    }

    return rotation0_found && rotation1_found;
}

void test_rotated_verticies(Sketch &s, std::vector<size_t> index, std::shared_ptr<Vertex> center, double angle, size_t copies) {
    for (size_t j = 1; j != (copies + 1); ++j) {
        for (size_t i : index) {
            EXPECT_TRUE(has_rotational_image(s, s.vertex(i), center, angle * j));
        }
    }
}

void test_rotated_curves(Sketch &s, std::vector<size_t> index, std::shared_ptr<Vertex> center, double angle, size_t copies) {
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