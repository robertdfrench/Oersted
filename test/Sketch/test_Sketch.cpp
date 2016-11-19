#include "test_Sketch.hpp"

void test_sketch_size(Sketch &s, size_t nverts, size_t ncurves, size_t nconstraints, size_t ncontours) {
    EXPECT_EQ(nverts, s.size_verticies());
    EXPECT_EQ(ncurves, s.size_curves());
    EXPECT_EQ(nconstraints, s.size_constraints());
    EXPECT_EQ(ncontours, s.size_contours());
}

TEST(Sketch, constructor) {
    Sketch s;
}

TEST(Sketch, new_element) {
    { // Case 0: Vertex
        Sketch s;

        auto v0 = s.new_element_SHARED_PTR<Vertex>();
        EXPECT_TRUE(s.size() == 1);
        EXPECT_TRUE(v0->get_equation_index() == 0);
        EXPECT_TRUE(v0->X->get_index() == 0);
        EXPECT_TRUE(v0->Y->get_index() == 1);

        auto v1 = s.new_element_SHARED_PTR<Vertex>();

        EXPECT_TRUE(s.size() == 2);
        EXPECT_TRUE(v1->get_equation_index() == 0);
        EXPECT_TRUE(v1->X->get_index() == 2);
        EXPECT_TRUE(v1->Y->get_index() == 3);
    }

    { // Case 1: LineSegment
        Sketch s;

        auto line = s.new_element_SHARED_PTR<LineSegment>();
        EXPECT_TRUE(s.size() == 1);
    }

    { // Case 2: CircularArc
        Sketch s;

        auto line = s.new_element_SHARED_PTR<LineSegment>();
        EXPECT_TRUE(s.size() == 1);
    }

    { // Case 3: Length
        Sketch s;

        auto v0 = s.new_element_SHARED_PTR<Vertex>(3.14159, 2.7183);
        auto v1 = s.new_element_SHARED_PTR<Vertex>(6.14159, 6.7183);

        auto line = s.new_element_SHARED_PTR<LineSegment>(v0, v1);

        auto length = s.new_element_SHARED_PTR<Length>(line, 1.0);

        EXPECT_TRUE(s.size() == 4);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Sketch__new_element_Length");
    }
}