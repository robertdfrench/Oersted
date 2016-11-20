#include "test_Sketch.hpp"

TEST(Curve, supremum) {
    Sketch s;

    auto origin = s.new_element<Vertex>(0.0, 0.0);
    auto v0 = s.new_element<Vertex>(1.0, 0.0);
    auto v1 = s.new_element<Vertex>(4.0, 0.0);
    auto v2 = s.new_element<Vertex>(M_SQRT2, M_SQRT2);
    auto v3 = s.new_element<Vertex>(M_SQRT1_2, M_SQRT1_2);

    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v3, v2);
    auto c0 = s.new_element<CircularArc>(v1, v2, origin, 4.0);
    auto c1 = s.new_element<CircularArc>(v0, v3, origin, 1.0);

    auto r0 = s.new_element<Radius>(c0, 4.0);
    auto r1 = s.new_element<Radius>(c1, 1.0);
    auto h0 = s.new_element<Horizontal>(l0);
    auto a0 = s.new_element<Angle>(l0, l1, 45.0);
    auto coin0 = s.new_element<Coincident<LineSegment>>(origin, l0);
    auto coin1 = s.new_element<Coincident<LineSegment>>(origin, l1);

    std::pair<double,double> sc0 = c0->supremum();
    std::pair<double,double> sl0 = l0->supremum();

    EXPECT_GT(sc0, sl0);

    s.solve();
    s.build();

    s.delete_me();
}