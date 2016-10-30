#include "test_Sketch.hpp"

TEST(Curve, supremum) {
    Sketch s;

    Vertex &origin = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v0 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(4.0, 0.0);
    Vertex &v2 = s.new_element<Vertex>(M_SQRT2, M_SQRT2);
    Vertex &v3 = s.new_element<Vertex>(M_SQRT1_2, M_SQRT1_2);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v3, v2);
    CircularArc &c0 = s.new_element<CircularArc>(v1, v2, origin, 4.0);
    CircularArc &c1 = s.new_element<CircularArc>(v0, v3, origin, 1.0);

    Radius &r0 = s.new_element<Radius>(c0, 4.0);
    Radius &r1 = s.new_element<Radius>(c1, 1.0);
    Horizontal &h0 = s.new_element<Horizontal>(l0);
    Angle &a0 = s.new_element<Angle>(l0, l1, 45.0);
    Coincident<LineSegment> &coin0 = s.new_element<Coincident<LineSegment>>(origin, l0);
    Coincident<LineSegment> &coin1 = s.new_element<Coincident<LineSegment>>(origin, l1);

    std::pair<double,double> sc0 = c0.supremum();
    std::pair<double,double> sl0 = l0.supremum();

    EXPECT_GT(sc0, sl0);

    s.solve();
    s.build();
}