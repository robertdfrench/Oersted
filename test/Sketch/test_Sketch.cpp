#include "test_Sketch.hpp"

TEST(Sketch, constructor) {
    Sketch s;
}

TEST(Sketch, new_element) {
    { // Case 0: Vertex
        Sketch s;

        Vertex &v0 = s.new_element<Vertex>();
        EXPECT_TRUE(s.size() == 1);
        EXPECT_TRUE(v0.get_equation_index() == 0);
        EXPECT_TRUE(v0.X->get_index() == 0);
        EXPECT_TRUE(v0.Y->get_index() == 1);

        Vertex &v1 = s.new_element<Vertex>();

        EXPECT_TRUE(s.size() == 2);
        EXPECT_TRUE(v1.get_equation_index() == 0);
        EXPECT_TRUE(v1.X->get_index() == 2);
        EXPECT_TRUE(v1.Y->get_index() == 3);
    }

    { // Case 1: LineSegment
        Sketch s;

        LineSegment &line = s.new_element<LineSegment>();
        EXPECT_TRUE(s.size() == 1);
    }

    { // Case 2: CircularArc
        Sketch s;

        LineSegment &line = s.new_element<LineSegment>();
        EXPECT_TRUE(s.size() == 1);
    }

    { // Case 3: Length
        Sketch s;

        Vertex &v0 = s.new_element<Vertex>(3.14159, 2.7183);
        Vertex &v1 = s.new_element<Vertex>(6.14159, 6.7183);

        LineSegment &line = s.new_element<LineSegment>(v0, v1);

        Length &length = s.new_element<Length>(line, 1.0);

        EXPECT_TRUE(s.size() == 4);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Sketch__new_element_Length");
    }
}