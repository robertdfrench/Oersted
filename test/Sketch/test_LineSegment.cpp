#include "test_Sketch.hpp"

TEST(LineSegment, constructor) {
    {   // ARGS::()
        EXPECT_NO_THROW(LineSegment l);
    }

    {   // ARGS::(Vertex,Vertex)
        std::shared_ptr<Vertex> v0, v1;
        EXPECT_NO_THROW(LineSegment l(v0, v1));
    }
}

TEST(LineSegment, length) {
    std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(3.14159, 2.7183);
    std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(6.14159, 6.7183);
    LineSegment line{v0, v1};
    EXPECT_NEAR(5.0, line.length(), TOL);
}

TEST(LineSegment, point) {
    std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(-1.0, 1.0);
    std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(2.0, -3.0);

    LineSegment l{v0, v1};

    for (size_t i = 0; i < 11; ++i) {
        double s = i / 10.0;

        double2 v = l.point(s);
        EXPECT_NEAR(v0->x() * (1.0 - s) + v1->x() * s, v.X, TOL);
        EXPECT_NEAR(v0->y() * (1.0 - s) + v1->y() * s, v.Y, TOL);
    }
}

TEST(LineSegment, tangent) {
    std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(-1.0, 1.0);
    std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(2.0, -3.0);

    LineSegment l{v0, v1};

    for (size_t i = 0; i < 11; ++i) {
        double s = i / 10.0;

        double2 v = l.tangent(s, true);
        EXPECT_NEAR(3.0 / 5.0, v.X, TOL);
        EXPECT_NEAR(-4.0 / 5.0, v.Y, TOL);

        v = l.tangent(s, false);
        EXPECT_NEAR(-3.0 / 5.0, v.X, TOL);
        EXPECT_NEAR(4.0 / 5.0, v.Y, TOL);
    }
}

TEST(LineSegment, a) {
    std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(-1.0, 1.0);
    std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(2.0, -3.0);

    LineSegment l{v0, v1};

    for (size_t i = 0; i < 11; ++i) {
        double s = i / 10.0;

        double ta = l.a(s, true);
        EXPECT_NEAR(atan2(-4.0, 3.0), ta, TOL);

        ta = l.a(s, false);
        EXPECT_NEAR(atan2(4.0, -3.0), ta, TOL);
    }
}

TEST(LineSegment, da) {
    std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(-1.0, 1.0);
    std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(2.0, -3.0);

    LineSegment l{v0, v1};

    for (size_t i = 0; i < 11; ++i) {
        double s = i / 10.0;

        double dta = l.da(s, true);
        EXPECT_NEAR(0.0, dta, TOL);

        dta = l.da(s, false);
        EXPECT_NEAR(0.0, dta, TOL);
    }
}

TEST(LineSegment, on_manifold) {
    {   //ARGS::(Vertex)
        Sketch s;

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(1.0, 1.0);
        auto v2 = s.new_element<Vertex>(1.1, 1.1);
        auto v3 = s.new_element<Vertex>(M_SQRT2, 0.0);
        auto v4 = s.new_element<Vertex>(0.0, M_SQRT2);
        auto v5 = s.new_element<Vertex>(0.5, 0.5);
        auto v6 = s.new_element<Vertex>(0.5, sqrt(2.0 - 0.25));

        auto l0 = s.new_element<LineSegment>(v0, v1);

        EXPECT_TRUE(l0->on_manifold(v0));
        EXPECT_TRUE(l0->on_manifold(v1));
        EXPECT_TRUE(l0->on_manifold(v2));
        EXPECT_FALSE(l0->on_manifold(v3));
        EXPECT_FALSE(l0->on_manifold(v4));
        EXPECT_TRUE(l0->on_manifold(v5));
        EXPECT_FALSE(l0->on_manifold(v6));
    }

    {   //ARGS::(Vertex,Vertex,double)
        Sketch s;

        auto vl0 = s.new_element<Vertex>(1.0, 1.0);
        auto vl1 = s.new_element<Vertex>(2.0, 2.0);

        auto l = s.new_element<LineSegment>(vl0, vl1);

        auto origin = s.new_element<Vertex>(0.5, 0.5);
        auto v0 = s.new_element<Vertex>(0.5 + M_SQRT1_2, 0.5);
        auto v1 = s.new_element<Vertex>(1.0, 0.5);
        auto v2 = s.new_element<Vertex>(2.0, 0.5);
        auto v3 = s.new_element<Vertex>(0.5 + 3.0 * M_SQRT1_2, 0.5);
        auto v4 = s.new_element<Vertex>(0.5 + 3.0, 0.5);

        double a = 44.0;
        EXPECT_FALSE(l->on_manifold(v0, origin, a));
        EXPECT_FALSE(l->on_manifold(v1, origin, a));
        EXPECT_FALSE(l->on_manifold(v2, origin, a));
        EXPECT_FALSE(l->on_manifold(v3, origin, a));
        EXPECT_FALSE(l->on_manifold(v4, origin, a));

        a += 1.0;
        EXPECT_TRUE(l->on_manifold(v0, origin, a));
        EXPECT_TRUE(l->on_manifold(v1, origin, a));
        EXPECT_TRUE(l->on_manifold(v2, origin, a));
        EXPECT_TRUE(l->on_manifold(v3, origin, a));
        EXPECT_TRUE(l->on_manifold(v4, origin, a));

        a += 1.0;
        EXPECT_FALSE(l->on_manifold(v0, origin, a));
        EXPECT_FALSE(l->on_manifold(v1, origin, a));
        EXPECT_FALSE(l->on_manifold(v2, origin, a));
        EXPECT_FALSE(l->on_manifold(v3, origin, a));
        EXPECT_FALSE(l->on_manifold(v4, origin, a));
    }
}

TEST(LineSegment, on_segment) {
    {   // ARGS::(Vertex)
        Sketch s;

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(1.0, 1.0);
        auto l = s.new_element<LineSegment>(v0, v1);

        auto von0 = s.new_element<Vertex>(0.3, 0.3);
        auto von1 = s.new_element<Vertex>(0.6, 0.6);

        auto voff0 = s.new_element<Vertex>(-0.3, -0.3);
        auto voff1 = s.new_element<Vertex>(1.3, 1.3);
        auto voff2 = s.new_element<Vertex>(0.6, 0.3);
        auto voff3 = s.new_element<Vertex>(-0.6, 0.3);

        EXPECT_TRUE(l->on_segment(von0));
        EXPECT_TRUE(l->on_segment(von1));

        EXPECT_FALSE(l->on_segment(voff0));
        EXPECT_FALSE(l->on_segment(voff1));
        EXPECT_FALSE(l->on_segment(voff2));
        EXPECT_FALSE(l->on_segment(voff3));
    }

    { // ARGS::(Vertex,Vertex,double)
        Sketch s;

        double angle = 45.0;

        auto origin = s.new_element<Vertex>(1.0, 1.0);

        auto vl0 = s.new_element<Vertex>(2.0, 2.0);
        auto vl1 = s.new_element<Vertex>(3.0, 3.0);

        auto l0 = s.new_element<LineSegment>(vl0, vl1);

        auto vt0 = s.new_element<Vertex>(1.0 + M_SQRT2, 1.0);
        auto vt1 = s.new_element<Vertex>(1.0 + 1.5 * M_SQRT2, 1.0);
        auto vt2 = s.new_element<Vertex>(1.0 + 2.0 * M_SQRT2, 1.0);

        auto vf0 = s.new_element<Vertex>(1.0 + 0.5 * M_SQRT2, 1.0);
        auto vf1 = s.new_element<Vertex>(1.0 + 2.5 * M_SQRT2, 1.0);

        EXPECT_TRUE(l0->on_segment(vt0, origin, angle));
        EXPECT_TRUE(l0->on_segment(vt1, origin, angle));
        EXPECT_TRUE(l0->on_segment(vt2, origin, angle));

        EXPECT_FALSE(l0->on_segment(vf0, origin, angle));
        EXPECT_FALSE(l0->on_segment(vf1, origin, angle));
    }
}

TEST(LineSegment, is_identical) {
    {   // ARGS::(LineSegment)
        Sketch s;

        auto v00 = s.new_element<Vertex>(0.0, 0.0);
        auto v01 = s.new_element<Vertex>(1.0, 1.0);
        auto l = s.new_element<LineSegment>(v00, v01);
        auto lr = s.new_element<LineSegment>(v01, v00);

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(1.0, 1.0);
        auto v2 = s.new_element<Vertex>(1.0, 0.0);
        auto v3 = s.new_element<Vertex>(0.0, 1.0);
        auto v4 = s.new_element<Vertex>(0.5, 0.5);

        auto l0 = s.new_element<LineSegment>(v0, v1);
        auto l1 = s.new_element<LineSegment>(v1, v0);

        EXPECT_TRUE(l->is_identical(l) == MatchOrientation::Forward);
        EXPECT_TRUE(l->is_identical(lr) == MatchOrientation::Reverse);
        EXPECT_TRUE(l->is_identical(l0) == MatchOrientation::Forward);
        EXPECT_TRUE(l->is_identical(l1) == MatchOrientation::Reverse);

        auto l2 = s.new_element<LineSegment>(v0, v2);
        auto l3 = s.new_element<LineSegment>(v0, v3);
        auto l4 = s.new_element<LineSegment>(v1, v2);
        auto l5 = s.new_element<LineSegment>(v1, v3);
        auto l6 = s.new_element<LineSegment>(v0, v4);
        auto l7 = s.new_element<LineSegment>(v1, v4);

        EXPECT_TRUE(l->is_identical(l2) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l3) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l4) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l5) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l6) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l7) == MatchOrientation::None);
    }

    {   // ARGS::(LineSegment,Vertex,double)
        Sketch s;

        auto origin = s.new_element<Vertex>(1.0, 1.0);

        auto v00 = s.new_element<Vertex>(M_SQRT2 + 1.0, M_SQRT2 + 1.0);
        auto v01 = s.new_element<Vertex>(1.5 * M_SQRT2 + 1.0, 1.5 * M_SQRT2 + 1.0);
        auto l = s.new_element<LineSegment>(v00, v01);

        auto v0 = s.new_element<Vertex>(2.0 + 1.0, 0.0 + 1.0);
        auto v1 = s.new_element<Vertex>(3.0 + 1.0, 0.0 + 1.0);
        auto v2 = s.new_element<Vertex>(0.0 + 1.0, 2.0 + 1.0);
        auto v3 = s.new_element<Vertex>(0.0 + 1.0, 3.0 + 1.0);

        auto l0 = s.new_element<LineSegment>(v0, v1);
        auto l1 = s.new_element<LineSegment>(v1, v0);
        auto l2 = s.new_element<LineSegment>(v2, v3);
        auto l3 = s.new_element<LineSegment>(v3, v2);

        EXPECT_TRUE(l->is_identical(l0, origin, +45.0) == MatchOrientation::Forward);
        EXPECT_TRUE(l->is_identical(l1, origin, +45.0) == MatchOrientation::Reverse);
        EXPECT_TRUE(l->is_identical(l2, origin, -45.0) == MatchOrientation::Forward);
        EXPECT_TRUE(l->is_identical(l3, origin, -45.0) == MatchOrientation::Reverse);

        EXPECT_TRUE(l->is_identical(l0, origin, -45.0) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l1, origin, -45.0) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l2, origin, +45.0) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l3, origin, +45.0) == MatchOrientation::None);

        auto v4 = s.new_element<Vertex>(2.5 + 1.0, 0.0 + 1.0);
        auto v5 = s.new_element<Vertex>(0.0 + 1.0, 2.5 + 1.0);

        auto l4 = s.new_element<LineSegment>(v0, v4);
        auto l5 = s.new_element<LineSegment>(v4, v0);
        auto l6 = s.new_element<LineSegment>(v1, v4);
        auto l7 = s.new_element<LineSegment>(v4, v1);
        auto l8 = s.new_element<LineSegment>(v2, v5);
        auto l9 = s.new_element<LineSegment>(v5, v2);
        auto l10 = s.new_element<LineSegment>(v3, v5);
        auto l11 = s.new_element<LineSegment>(v5, v3);

        EXPECT_TRUE(l->is_identical(l4, origin, +45.0) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l5, origin, +45.0) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l6, origin, +45.0) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l7, origin, +45.0) == MatchOrientation::None);

        EXPECT_TRUE(l->is_identical(l8, origin, -45.0) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l9, origin, -45.0) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l10, origin, -45.0) == MatchOrientation::None);
        EXPECT_TRUE(l->is_identical(l11, origin, -45.0) == MatchOrientation::None);
    }
}

TEST(LineSegment, s_coincident) {
    {   // ARGS::(LineSegment)
        Sketch s;

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(1.0, 1.0);
        auto v2 = s.new_element<Vertex>(1.1, 1.1);
        auto v3 = s.new_element<Vertex>(2.0, 2.0);
        auto v4 = s.new_element<Vertex>(3.1, 0.9);
        auto v5 = s.new_element<Vertex>(M_SQRT2, 0.0);
        auto v6 = s.new_element<Vertex>(0.0, M_SQRT2);
        auto v7 = s.new_element<Vertex>(2.0 * M_SQRT2, 0.0);
        auto v8 = s.new_element<Vertex>(1.0, 1.0 - M_SQRT2);
        auto v9 = s.new_element<Vertex>(M_SQRT2 + 1.0, 1.0);

        auto l0 = s.new_element<LineSegment>(v0, v1);
        auto l1 = s.new_element<LineSegment>(v2, v3);
        auto l2 = s.new_element<LineSegment>(v3, v2);
        auto l3 = s.new_element<LineSegment>(v2, v4);
        auto l4 = s.new_element<LineSegment>(v4, v2);
        auto l5 = s.new_element<LineSegment>(v5, v9);

        EXPECT_TRUE(l0->is_coincident(l0));
        EXPECT_TRUE(l0->is_coincident(l1));
        EXPECT_TRUE(l0->is_coincident(l2));
        EXPECT_FALSE(l0->is_coincident(l3));
        EXPECT_FALSE(l0->is_coincident(l4));
        EXPECT_FALSE(l0->is_coincident(l5));
    }

    {   // ARGS::(CircularArc)
        LineSegment l = LineSegment();
        std::shared_ptr<CircularArc> c = std::make_shared<CircularArc>();

        EXPECT_FALSE(l.is_coincident(c));
    }
}