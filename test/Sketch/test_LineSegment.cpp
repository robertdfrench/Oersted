#include "test_Sketch.hpp"

TEST(LINESEGMENT, CONSTRUCTOR) {
    {   // ARGS::()
        EXPECT_NO_THROW(LineSegment l);
    }

    {   // ARGS::(Vertex,Vertex)
        Vertex v0, v1;
        EXPECT_NO_THROW(LineSegment l(v0, v1));
    }
}

TEST(LINESEGMENT, METHOD_length) {
    Vertex v0{3.14159, 2.7183};
    Vertex v1{6.14159, 6.7183};
    LineSegment line{v0, v1};
    EXPECT_NEAR(5.0, line.length(), TOL);
}

TEST(LINESEGMENT, METHOD_point) {
    Vertex v0{-1.0, 1.0};
    Vertex v1{2.0, -3.0};

    LineSegment l{v0, v1};

    for (size_t i = 0; i < 11; ++i) {
        double s = i / 10.0;

        Vertex v = l.point(s);
        EXPECT_NEAR(v0.x() * (1.0 - s) + v1.x() * s, v.x(), TOL);
        EXPECT_NEAR(v0.y() * (1.0 - s) + v1.y() * s, v.y(), TOL);
    }
}

TEST(LINESEGMENT, METHOD_tangent) {
    Vertex v0{-1.0, 1.0};
    Vertex v1{2.0, -3.0};

    LineSegment l{v0, v1};

    for (size_t i = 0; i < 11; ++i) {
        double s = i / 10.0;

        Vertex v = l.tangent(s, true);
        EXPECT_NEAR(3.0 / 5.0, v.x(), TOL);
        EXPECT_NEAR(-4.0 / 5.0, v.y(), TOL);

        v = l.tangent(s, false);
        EXPECT_NEAR(-3.0 / 5.0, v.x(), TOL);
        EXPECT_NEAR(4.0 / 5.0, v.y(), TOL);
    }
}

TEST(LINESEGMENT, METHOD_a) {
    Vertex v0{-1.0, 1.0};
    Vertex v1{2.0, -3.0};

    LineSegment l{v0, v1};

    for (size_t i = 0; i < 11; ++i) {
        double s = i / 10.0;

        double ta = l.a(s, true);
        EXPECT_NEAR(atan2(-4.0, 3.0), ta, TOL);

        ta = l.a(s, false);
        EXPECT_NEAR(atan2(4.0, -3.0), ta, TOL);
    }
}

TEST(LINESEGMENT, METHOD_da) {
    Vertex v0{-1.0, 1.0};
    Vertex v1{2.0, -3.0};

    LineSegment l{v0, v1};

    for (size_t i = 0; i < 11; ++i) {
        double s = i / 10.0;

        double dta = l.da(s, true);
        EXPECT_NEAR(0.0, dta, TOL);

        dta = l.da(s, false);
        EXPECT_NEAR(0.0, dta, TOL);
    }
}

TEST(LINESEGMENT, METHOD_on_manifold) {
    {   //ARGS::(Vertex)
        Sketch s;

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(1.0, 1.0);
        Vertex &v2 = s.new_element<Vertex>(1.1, 1.1);
        Vertex &v3 = s.new_element<Vertex>(M_SQRT2, 0.0);
        Vertex &v4 = s.new_element<Vertex>(0.0, M_SQRT2);
        Vertex &v5 = s.new_element<Vertex>(0.5, 0.5);
        Vertex &v6 = s.new_element<Vertex>(0.5, sqrt(2.0 - 0.25));

        LineSegment &l0 = s.new_element<LineSegment>(v0, v1);

        EXPECT_TRUE(l0.on_manifold(&v0));
        EXPECT_TRUE(l0.on_manifold(&v1));
        EXPECT_TRUE(l0.on_manifold(&v2));
        EXPECT_FALSE(l0.on_manifold(&v3));
        EXPECT_FALSE(l0.on_manifold(&v4));
        EXPECT_TRUE(l0.on_manifold(&v5));
        EXPECT_FALSE(l0.on_manifold(&v6));
    }

    {   //ARGS::(Vertex,Vertex,double)
        Sketch s;

        Vertex &vl0 = s.new_element<Vertex>(1.0, 1.0);
        Vertex &vl1 = s.new_element<Vertex>(2.0, 2.0);

        LineSegment &l = s.new_element<LineSegment>(vl0, vl1);

        Vertex &origin = s.new_element<Vertex>(0.5, 0.5);
        Vertex &v0 = s.new_element<Vertex>(0.5 + M_SQRT1_2, 0.5);
        Vertex &v1 = s.new_element<Vertex>(1.0, 0.5);
        Vertex &v2 = s.new_element<Vertex>(2.0, 0.5);
        Vertex &v3 = s.new_element<Vertex>(0.5 + 3.0 * M_SQRT1_2, 0.5);
        Vertex &v4 = s.new_element<Vertex>(0.5 + 3.0, 0.5);

        double a = 44.0;
        EXPECT_FALSE(l.on_manifold(&v0, &origin, a));
        EXPECT_FALSE(l.on_manifold(&v1, &origin, a));
        EXPECT_FALSE(l.on_manifold(&v2, &origin, a));
        EXPECT_FALSE(l.on_manifold(&v3, &origin, a));
        EXPECT_FALSE(l.on_manifold(&v4, &origin, a));

        a += 1.0;
        EXPECT_TRUE(l.on_manifold(&v0, &origin, a));
        EXPECT_TRUE(l.on_manifold(&v1, &origin, a));
        EXPECT_TRUE(l.on_manifold(&v2, &origin, a));
        EXPECT_TRUE(l.on_manifold(&v3, &origin, a));
        EXPECT_TRUE(l.on_manifold(&v4, &origin, a));

        a += 1.0;
        EXPECT_FALSE(l.on_manifold(&v0, &origin, a));
        EXPECT_FALSE(l.on_manifold(&v1, &origin, a));
        EXPECT_FALSE(l.on_manifold(&v2, &origin, a));
        EXPECT_FALSE(l.on_manifold(&v3, &origin, a));
        EXPECT_FALSE(l.on_manifold(&v4, &origin, a));
    }
}

TEST(LINESEGMENT, METHOD_on_segment) {
    {   // ARGS::(Vertex)
        Sketch s;

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(1.0, 1.0);
        LineSegment &l = s.new_element<LineSegment>(v0, v1);

        Vertex &von0 = s.new_element<Vertex>(0.3, 0.3);
        Vertex &von1 = s.new_element<Vertex>(0.6, 0.6);

        Vertex &voff0 = s.new_element<Vertex>(-0.3, -0.3);
        Vertex &voff1 = s.new_element<Vertex>(1.3, 1.3);
        Vertex &voff2 = s.new_element<Vertex>(0.6, 0.3);
        Vertex &voff3 = s.new_element<Vertex>(-0.6, 0.3);

        EXPECT_TRUE(l.on_segment(&von0));
        EXPECT_TRUE(l.on_segment(&von1));

        EXPECT_FALSE(l.on_segment(&voff0));
        EXPECT_FALSE(l.on_segment(&voff1));
        EXPECT_FALSE(l.on_segment(&voff2));
        EXPECT_FALSE(l.on_segment(&voff3));
    }

    { // ARGS::(Vertex,Vertex,double)
        Sketch s;

        double angle = 45.0;

        Vertex &origin = s.new_element<Vertex>(1.0, 1.0);

        Vertex &vl0 = s.new_element<Vertex>(2.0, 2.0);
        Vertex &vl1 = s.new_element<Vertex>(3.0, 3.0);

        LineSegment &l0 = s.new_element<LineSegment>(vl0, vl1);

        Vertex &vt0 = s.new_element<Vertex>(1.0 + M_SQRT2, 1.0);
        Vertex &vt1 = s.new_element<Vertex>(1.0 + 1.5 * M_SQRT2, 1.0);
        Vertex &vt2 = s.new_element<Vertex>(1.0 + 2.0 * M_SQRT2, 1.0);

        Vertex &vf0 = s.new_element<Vertex>(1.0 + 0.5 * M_SQRT2, 1.0);
        Vertex &vf1 = s.new_element<Vertex>(1.0 + 2.5 * M_SQRT2, 1.0);

        EXPECT_TRUE(l0.on_segment(&vt0, &origin, angle));
        EXPECT_TRUE(l0.on_segment(&vt1, &origin, angle));
        EXPECT_TRUE(l0.on_segment(&vt2, &origin, angle));

        EXPECT_FALSE(l0.on_segment(&vf0, &origin, angle));
        EXPECT_FALSE(l0.on_segment(&vf1, &origin, angle));
    }
}

TEST(LINESEGMENT, METHOD_is_identical) {
    {   // ARGS::(LineSegment)
        Sketch s;

        Vertex &v00 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v01 = s.new_element<Vertex>(1.0, 1.0);
        LineSegment &l = s.new_element<LineSegment>(v00, v01);
        LineSegment &lr = s.new_element<LineSegment>(v01, v00);

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(1.0, 1.0);
        Vertex &v2 = s.new_element<Vertex>(1.0, 0.0);
        Vertex &v3 = s.new_element<Vertex>(0.0, 1.0);
        Vertex &v4 = s.new_element<Vertex>(0.5, 0.5);

        LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
        LineSegment &l1 = s.new_element<LineSegment>(v1, v0);

        EXPECT_TRUE(l.is_identical(&l));
        EXPECT_TRUE(l.is_identical(&lr));
        EXPECT_TRUE(l.is_identical(&l0));
        EXPECT_TRUE(l.is_identical(&l1));

        LineSegment &l2 = s.new_element<LineSegment>(v0, v2);
        LineSegment &l3 = s.new_element<LineSegment>(v0, v3);
        LineSegment &l4 = s.new_element<LineSegment>(v1, v2);
        LineSegment &l5 = s.new_element<LineSegment>(v1, v3);
        LineSegment &l6 = s.new_element<LineSegment>(v0, v4);
        LineSegment &l7 = s.new_element<LineSegment>(v1, v4);

        EXPECT_FALSE(l.is_identical(&l2));
        EXPECT_FALSE(l.is_identical(&l3));
        EXPECT_FALSE(l.is_identical(&l4));
        EXPECT_FALSE(l.is_identical(&l5));
        EXPECT_FALSE(l.is_identical(&l6));
        EXPECT_FALSE(l.is_identical(&l7));
    }

    {   // ARGS::(LineSegment,Vertex,double)
        Sketch s;

        Vertex &origin = s.new_element<Vertex>(1.0, 1.0);

        Vertex &v00 = s.new_element<Vertex>(M_SQRT2 + 1.0, M_SQRT2 + 1.0);
        Vertex &v01 = s.new_element<Vertex>(1.5 * M_SQRT2 + 1.0, 1.5 * M_SQRT2 + 1.0);
        LineSegment &l = s.new_element<LineSegment>(v00, v01);

        Vertex &v0 = s.new_element<Vertex>(2.0 + 1.0, 0.0 + 1.0);
        Vertex &v1 = s.new_element<Vertex>(3.0 + 1.0, 0.0 + 1.0);
        Vertex &v2 = s.new_element<Vertex>(0.0 + 1.0, 2.0 + 1.0);
        Vertex &v3 = s.new_element<Vertex>(0.0 + 1.0, 3.0 + 1.0);

        LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
        LineSegment &l1 = s.new_element<LineSegment>(v1, v0);
        LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
        LineSegment &l3 = s.new_element<LineSegment>(v3, v2);

        EXPECT_TRUE(l.is_identical(&l0, &origin, +45.0));
        EXPECT_TRUE(l.is_identical(&l1, &origin, +45.0));
        EXPECT_TRUE(l.is_identical(&l2, &origin, -45.0));
        EXPECT_TRUE(l.is_identical(&l3, &origin, -45.0));

        EXPECT_FALSE(l.is_identical(&l0, &origin, -45.0));
        EXPECT_FALSE(l.is_identical(&l1, &origin, -45.0));
        EXPECT_FALSE(l.is_identical(&l2, &origin, +45.0));
        EXPECT_FALSE(l.is_identical(&l3, &origin, +45.0));

        Vertex &v4 = s.new_element<Vertex>(2.5 + 1.0, 0.0 + 1.0);
        Vertex &v5 = s.new_element<Vertex>(0.0 + 1.0, 2.5 + 1.0);

        LineSegment &l4 = s.new_element<LineSegment>(v0, v4);
        LineSegment &l5 = s.new_element<LineSegment>(v4, v0);
        LineSegment &l6 = s.new_element<LineSegment>(v1, v4);
        LineSegment &l7 = s.new_element<LineSegment>(v4, v1);
        LineSegment &l8 = s.new_element<LineSegment>(v2, v5);
        LineSegment &l9 = s.new_element<LineSegment>(v5, v2);
        LineSegment &l10 = s.new_element<LineSegment>(v3, v5);
        LineSegment &l11 = s.new_element<LineSegment>(v5, v3);

        EXPECT_FALSE(l.is_identical(&l4, &origin, +45.0));
        EXPECT_FALSE(l.is_identical(&l5, &origin, +45.0));
        EXPECT_FALSE(l.is_identical(&l6, &origin, +45.0));
        EXPECT_FALSE(l.is_identical(&l7, &origin, +45.0));

        EXPECT_FALSE(l.is_identical(&l8, &origin, -45.0));
        EXPECT_FALSE(l.is_identical(&l9, &origin, -45.0));
        EXPECT_FALSE(l.is_identical(&l10, &origin, -45.0));
        EXPECT_FALSE(l.is_identical(&l11, &origin, -45.0));
    }
}

TEST(LINESEGMENT, METHOD_is_coincident) {
    {   // ARGS::(LineSegment)
        Sketch s;

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(1.0, 1.0);
        Vertex &v2 = s.new_element<Vertex>(1.1, 1.1);
        Vertex &v3 = s.new_element<Vertex>(2.0, 2.0);
        Vertex &v4 = s.new_element<Vertex>(3.1, 0.9);
        Vertex &v5 = s.new_element<Vertex>(M_SQRT2, 0.0);
        Vertex &v6 = s.new_element<Vertex>(0.0, M_SQRT2);
        Vertex &v7 = s.new_element<Vertex>(2.0 * M_SQRT2, 0.0);
        Vertex &v8 = s.new_element<Vertex>(1.0, 1.0 - M_SQRT2);
        Vertex &v9 = s.new_element<Vertex>(M_SQRT2 + 1.0, 1.0);

        LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
        LineSegment &l1 = s.new_element<LineSegment>(v2, v3);
        LineSegment &l2 = s.new_element<LineSegment>(v3, v2);
        LineSegment &l3 = s.new_element<LineSegment>(v2, v4);
        LineSegment &l4 = s.new_element<LineSegment>(v4, v2);
        LineSegment &l5 = s.new_element<LineSegment>(v5, v9);

        EXPECT_TRUE(l0.is_coincident(&l0));
        EXPECT_TRUE(l0.is_coincident(&l1));
        EXPECT_TRUE(l0.is_coincident(&l2));
        EXPECT_FALSE(l0.is_coincident(&l3));
        EXPECT_FALSE(l0.is_coincident(&l4));
        EXPECT_FALSE(l0.is_coincident(&l5));
    }

    {   // ARGS::(CircularArc)
        LineSegment l = LineSegment();
        CircularArc c = CircularArc();

        EXPECT_FALSE(l.is_coincident(&c));
    }
}