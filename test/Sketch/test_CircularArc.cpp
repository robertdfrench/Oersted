#include <math.h>
#include <tgmath.h>
#include "test_Sketch.hpp"

TEST(CIRCULARARC, CONSTRUCTOR) {
    { //ARGS::()
        EXPECT_NO_THROW(CircularArc c);
    }

    { //ARGS::(Vertex,Vertex,Vertex)
        Vertex v0, v1, center;
        EXPECT_NO_THROW(CircularArc c(v0, v1, center));
    }
}

TEST(CIRCULARARC, METHOD_point) {
    {
        double r = 1.0;
        double a0 = -M_PI_2;
        double a1 = M_PI_2;

        Vertex vc{0.0, 0.0};
        Vertex v0{r * cos(a0), r * sin(a0)};
        Vertex v1{r * cos(a1), r * sin(a1)};
        CircularArc c{v0, v1, vc, r};

        for (size_t i = 0; i < 11; ++i) {
            double s = i / 10.0;
            double a = a0 * (1.0 - s) + a1 * s;

            Vertex v = c.point(s);
            EXPECT_NEAR(r * cos(a), v.x(), TOL);
            EXPECT_NEAR(r * sin(a), v.y(), TOL);
        }
    }

    {
        double r = 1.0;
        double a0 = M_PI_2;
        double a1 = M_PI_2 + M_PI;

        Vertex vc{0.0, 0.0};
        Vertex v0{r * cos(a0), r * sin(a0)};
        Vertex v1{r * cos(a1), r * sin(a1)};
        CircularArc c{v0, v1, vc, r};

        for (size_t i = 0; i < 11; ++i) {
            double s = i / 10.0;
            double a = a0 * (1.0 - s) + a1 * s;

            Vertex v = c.point(s);
            EXPECT_NEAR(r * cos(a), v.x(), TOL);
            EXPECT_NEAR(r * sin(a), v.y(), TOL);
        }
    }
}

TEST(CIRCULARARC, METHOD_tangent) {
    {
        double r = 1.0;
        double a0 = -M_PI_2;
        double a1 = M_PI_2;

        Vertex vc{0.0, 0.0};
        Vertex v0{r * cos(a0), r * sin(a0)};
        Vertex v1{r * cos(a1), r * sin(a1)};
        CircularArc c{v0, v1, vc, r};

        for (size_t i = 0; i < 11; ++i) {
            double s = i / 10.0;
            double a = a0 * (1.0 - s) + a1 * s;

            Vertex v = c.tangent(s, true);
            EXPECT_NEAR(-sin(a), v.x(), TOL);
            EXPECT_NEAR(cos(a), v.y(), TOL);

            v = c.tangent(s, false);
            EXPECT_NEAR(sin(a), v.x(), TOL);
            EXPECT_NEAR(-cos(a), v.y(), TOL);
        }
    }

    {
        double r = 1.0;
        double a0 = M_PI_2;
        double a1 = M_PI_2 + M_PI;

        Vertex vc{0.0, 0.0};
        Vertex v0{r * cos(a0), r * sin(a0)};
        Vertex v1{r * cos(a1), r * sin(a1)};
        CircularArc c{v0, v1, vc, r};

        for (size_t i = 0; i < 11; ++i) {
            double s = i / 10.0;
            double a = a0 * (1.0 - s) + a1 * s;

            Vertex v = c.tangent(s, true);
            EXPECT_NEAR(-sin(a), v.x(), TOL);
            EXPECT_NEAR(cos(a), v.y(), TOL);

            v = c.tangent(s, false);
            EXPECT_NEAR(sin(a), v.x(), TOL);
            EXPECT_NEAR(-cos(a), v.y(), TOL);
        }
    }
}

TEST(CIRCULARARC, METHOD_a) {
    {
        double r = 1.0;
        double a0 = -M_PI_2;
        double a1 = M_PI_2;

        Vertex vc{0.0, 0.0};
        Vertex v0{r * cos(a0), r * sin(a0)};
        Vertex v1{r * cos(a1), r * sin(a1)};
        CircularArc c{v0, v1, vc, r};

        for (size_t i = 0; i < 11; ++i) {
            double s = i / 10.0;
            double a = a0 * (1.0 - s) + a1 * s;

            double ta = c.a(s, true);
            if (i == 10) {
                // Tangent derivative hits branch cut, ta == +- M_PI
                EXPECT_NEAR(M_PI, abs(ta), TOL);
            } else {
                EXPECT_NEAR(atan2(cos(a), -sin(a)), ta, TOL);
            }

            ta = c.a(s, false);
            if (i == 0) {
                // Tangent derivative hits branch cut, ta == +- M_PI
                EXPECT_NEAR(M_PI, abs(ta), TOL);
            } else {
                EXPECT_NEAR(atan2(-cos(a), sin(a)), ta, TOL);
            }
        }
    }

    {
        double r = 1.0;
        double a0 = M_PI_2;
        double a1 = M_PI_2 + M_PI;

        Vertex vc{0.0, 0.0};
        Vertex v0{r * cos(a0), r * sin(a0)};
        Vertex v1{r * cos(a1), r * sin(a1)};
        CircularArc c{v0, v1, vc, r};

        for (size_t i = 0; i < 11; ++i) {
            double s = i / 10.0;
            double a = a0 * (1.0 - s) + a1 * s;

            double ta = c.a(s, true);
            if (i == 0) {
                // Tangent derivative hits branch cut, ta == +- M_PI
                EXPECT_NEAR(M_PI, abs(ta), TOL);
            } else {
                EXPECT_NEAR(atan2(cos(a), -sin(a)), ta, TOL);
            }

            ta = c.a(s, false);
            if (i == 10) {
                // Tangent derivative hits branch cut, ta == +- M_PI
                EXPECT_NEAR(M_PI, abs(ta), TOL);
            } else {
                EXPECT_NEAR(atan2(-cos(a), sin(a)), ta, TOL);
            }
        }
    }
}

TEST(CIRCULARARC, METHOD_da) {
    {
        double r = 1.0;
        double a0 = -M_PI_2;
        double a1 = M_PI_2;

        Vertex vc{0.0, 0.0};
        Vertex v0{r * cos(a0), r * sin(a0)};
        Vertex v1{r * cos(a1), r * sin(a1)};
        CircularArc c{v0, v1, vc, r};

        for (size_t i = 0; i < 11; ++i) {
            double s = i / 10.0;
            double a = a0 * (1.0 - s) + a1 * s;

            double dta = c.da(s, true);
            EXPECT_NEAR(1.0 / r, dta, TOL);

            dta = c.da(s, false);
            EXPECT_NEAR(-1.0 / r, dta, TOL);
        }
    }

    {
        double r = 1.0;
        double a0 = M_PI_2;
        double a1 = M_PI_2 + M_PI;

        Vertex vc{0.0, 0.0};
        Vertex v0{r * cos(a0), r * sin(a0)};
        Vertex v1{r * cos(a1), r * sin(a1)};
        CircularArc c{v0, v1, vc, r};

        for (size_t i = 0; i < 11; ++i) {
            double s = i / 10.0;
            double a = a0 * (1.0 - s) + a1 * s;

            double dta = c.da(s, true);
            EXPECT_NEAR(1.0 / r, dta, TOL);

            dta = c.da(s, false);
            EXPECT_NEAR(-1.0 / r, dta, TOL);
        }
    }
}

TEST(CIRCULARARC, METHOD_on_manifold) {
    { //ARGS::(Vertex)
        Sketch s;

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(1.0, 1.0);
        Vertex &v2 = s.new_element<Vertex>(1.1, 1.1);
        Vertex &v3 = s.new_element<Vertex>(M_SQRT2, 0.0);
        Vertex &v4 = s.new_element<Vertex>(0.0, M_SQRT2);
        Vertex &v5 = s.new_element<Vertex>(0.5, 0.5);
        Vertex &v6 = s.new_element<Vertex>(0.5, sqrt(2.0 - 0.25));

        CircularArc &c0 = s.new_element<CircularArc>(v3, v1, v0, M_SQRT2);

        EXPECT_FALSE(c0.on_manifold(&v0));
        EXPECT_TRUE(c0.on_manifold(&v1));
        EXPECT_FALSE(c0.on_manifold(&v2));
        EXPECT_TRUE(c0.on_manifold(&v3));
        EXPECT_TRUE(c0.on_manifold(&v4));
        EXPECT_FALSE(c0.on_manifold(&v5));
        EXPECT_TRUE(c0.on_manifold(&v6));
    }

    { //ARGS::(Vertex,Vertex,double)
        Sketch s;

        Vertex &vcc = s.new_element<Vertex>(0.0, 0.0);
        Vertex &vc0 = s.new_element<Vertex>(1.0, 0.0);
        Vertex &vc1 = s.new_element<Vertex>(0.0, 1.0);
        double rc = 1.0;

        CircularArc &c = s.new_element<CircularArc>(vc0, vc1, vcc, rc);

        Vertex &v0 = s.new_element<Vertex>(M_SQRT1_2, M_SQRT1_2);
        double a0 = 90.0;
        EXPECT_TRUE(c.on_manifold(&v0, &vcc, a0));
        EXPECT_TRUE(c.on_manifold(&v0, &vcc, a0 + 1.0));
        EXPECT_TRUE(c.on_manifold(&v0, &vcc, a0 - 1.0));

        double r1 = hypot(M_SQRT1_2 - 1.0, M_SQRT1_2);
        double a1 = atan2(M_SQRT1_2, M_SQRT1_2 - 1.0) * 180.0 / M_PI;
        Vertex &v1 = s.new_element<Vertex>(1.0 + r1, 0.0);

        EXPECT_TRUE(c.on_manifold(&v1, &vc0, a1));
        EXPECT_FALSE(c.on_manifold(&v1, &vc0, a1 + 1.0));
        EXPECT_FALSE(c.on_manifold(&v1, &vc0, a1 - 1.0));

        Vertex &v2 = s.new_element<Vertex>(0.0, 3.0);
        Vertex &vo2 = s.new_element<Vertex>(0.0, 2.0);
        double a2 = -180.0;

        EXPECT_TRUE(c.on_manifold(&v2, &vo2, a2));
        EXPECT_FALSE(c.on_manifold(&v2, &vo2, a2 + 1.0));
        EXPECT_FALSE(c.on_manifold(&v2, &vo2, a2 - 1.0));

        Vertex &v3 = s.new_element<Vertex>(3.0, 0.0);
        Vertex &vo3 = s.new_element<Vertex>(2.0, 0.0);
        double a3 = 180.0;

        EXPECT_TRUE(c.on_manifold(&v3, &vo3, a3));
        EXPECT_FALSE(c.on_manifold(&v3, &vo3, a3 + 1.0));
        EXPECT_FALSE(c.on_manifold(&v3, &vo3, a3 - 1.0));
    }
}

TEST(CIRCULARARC, METHOD_on_segment) {
    { //ARGS::(Vertex)
        Sketch s;

        double r = 1.0;
        Vertex &v0 = s.new_element<Vertex>(1.0, 2.0);
        Vertex &v1 = s.new_element<Vertex>(1.0 + r, 2.0);
        Vertex &v2 = s.new_element<Vertex>(1.0, 2.0 + r);

        CircularArc &c = s.new_element<CircularArc>(v1, v2, v0, r);

        Vertex &von1 = v1;
        Vertex &von2 = v2;
        Vertex &von3 = s.new_element<Vertex>(1.0 + r * cos(M_PI_4), 2.0 + r * sin(M_PI_4));

        Vertex &voff1 = v0;
        Vertex &voff2 = s.new_element<Vertex>(1.0 + r * cos(-M_PI_4), 2.0 + r * sin(-M_PI_4));
        Vertex &voff3 = s.new_element<Vertex>(1.0 + r * cos(3.0 * M_PI_4), 2.0 + r * sin(3.0 * M_PI_4));
        Vertex &voff4 = s.new_element<Vertex>(1.0 + 2.0 * r * cos(M_PI_4), 2.0 + 2.0 * r * sin(M_PI_4));

        EXPECT_TRUE(c.on_segment(&von1));
        EXPECT_TRUE(c.on_segment(&von2));
        EXPECT_TRUE(c.on_segment(&von3));

        EXPECT_FALSE(c.on_segment(&voff1));
        EXPECT_FALSE(c.on_segment(&voff2));
        EXPECT_FALSE(c.on_segment(&voff3));
        EXPECT_FALSE(c.on_segment(&voff4));
    }

    {   // ARGS::(Vertex,Vertex,double)
        Sketch s;

        double r = 1.0;
        Vertex &v0 = s.new_element<Vertex>(1.0, 2.0);
        Vertex &v1 = s.new_element<Vertex>(1.0 + r, 2.0);
        Vertex &v2 = s.new_element<Vertex>(1.0, 2.0 + r);

        CircularArc &c = s.new_element<CircularArc>(v1, v2, v0, r);

        Vertex &origin = s.new_element<Vertex>(0.0, 0.0);

        Vertex &von1 = s.new_element<Vertex>(1.0 + r, -2.0);
        double aon1 = 2.0 * atan2(-von1.y(), von1.x()) * 180.0 / M_PI;

        Vertex &von2 = s.new_element<Vertex>(1.0, -2.0 - r);
        double aon2 = 2.0 * atan2(-von2.y(), von2.x()) * 180.0 / M_PI;

        Vertex &von3 = s.new_element<Vertex>(1.0 + r * cos(M_PI_4), -2.0 - r * sin(M_PI_4));
        double aon3 = 2.0 * atan2(-von3.y(), von3.x()) * 180.0 / M_PI;

        Vertex &voff1 = s.new_element<Vertex>(1.0, -2.0);
        double aoff1 = 2.0 * atan2(-voff1.y(), voff1.x()) * 180.0 / M_PI;

        Vertex &voff2 = s.new_element<Vertex>(1.0 + r * cos(-M_PI_4), -2.0 - r * sin(-M_PI_4));
        double aoff2 = 2.0 * atan2(-voff2.y(), voff2.x()) * 180.0 / M_PI;

        Vertex &voff3 = s.new_element<Vertex>(1.0 + r * cos(3.0 * M_PI_4), -2.0 - r * sin(3.0 * M_PI_4));
        double aoff3 = 2.0 * atan2(-voff3.y(), voff3.x()) * 180.0 / M_PI;

        Vertex &voff4 = s.new_element<Vertex>(1.0 + 2.0 * r * cos(M_PI_4), -2.0 - 2.0 * r * sin(M_PI_4));
        double aoff4 = 2.0 * atan2(-voff4.y(), voff4.x()) * 180.0 / M_PI;

        EXPECT_TRUE(c.on_segment(&von1, &origin, aon1));
        EXPECT_TRUE(c.on_segment(&von2, &origin, aon2));
        EXPECT_TRUE(c.on_segment(&von3, &origin, aon3));

        EXPECT_FALSE(c.on_segment(&voff1, &origin, aoff1));
        EXPECT_FALSE(c.on_segment(&voff2, &origin, aoff2));
        EXPECT_FALSE(c.on_segment(&voff3, &origin, aoff3));
        EXPECT_FALSE(c.on_segment(&voff4, &origin, aoff4));
    }
}

TEST(CIRCULARARC, METHOD_is_identical) {
    {   //ARGS::(Vertex)
        Sketch s;

        Vertex &vcc = s.new_element<Vertex>(0.0, 0.0);
        Vertex &vc0 = s.new_element<Vertex>(1.0, 0.0);
        Vertex &vc1 = s.new_element<Vertex>(0.0, 1.0);

        CircularArc &c = s.new_element<CircularArc>(vc0, vc1, vcc, 1.0);

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(1.0, 0.0);
        Vertex &v2 = s.new_element<Vertex>(0.0, 1.0);
        Vertex &v3 = s.new_element<Vertex>(1.0, 1.0);
        Vertex &v4 = s.new_element<Vertex>(M_SQRT1_2, M_SQRT1_2);

        CircularArc &c0 = s.new_element<CircularArc>(v1, v2, v0, 1.0);

        // True
        EXPECT_TRUE(c.is_identical(&c));
        EXPECT_TRUE(c.is_identical(&c0));

        // #TODO:	Radius does not match distance of endpoints from center
        //			Could be identical depending on other constraints
        //			Behavior is undefined unless Sketch::solve() is called
        CircularArc &c1 = s.new_element<CircularArc>(vc0, vc1, vcc, 0.5);
        EXPECT_FALSE(c.is_identical(&c1));

        CircularArc &c2 = s.new_element<CircularArc>(v1, v2, v0, 0.5);
        EXPECT_FALSE(c.is_identical(&c2));

        // False
        CircularArc &c3 = s.new_element<CircularArc>(vc1, vc0, vcc, 0.5);
        CircularArc &c4 = s.new_element<CircularArc>(v2, v1, v0, 1.0);
        CircularArc &c5 = s.new_element<CircularArc>(v1, v2, v3, 1.0);
        CircularArc &c6 = s.new_element<CircularArc>(v2, v1, v3, 1.0);

        EXPECT_FALSE(c.is_identical(&c3));
        EXPECT_FALSE(c.is_identical(&c4));
        EXPECT_FALSE(c.is_identical(&c5));
        EXPECT_FALSE(c.is_identical(&c6));
    }

    {   // ARGS::(Vertex,Vertex,double)
        Sketch s;

        Vertex &vc = s.new_element<Vertex>(1.0, 1.0);
        Vertex &vs = s.new_element<Vertex>(0.0, 1.0);
        Vertex &ve = s.new_element<Vertex>(1.0, 0.0);

        CircularArc &c = s.new_element<CircularArc>(vs, ve, vc, 1.0);

        Vertex &vc0 = s.new_element<Vertex>(1.0, 1.0);
        Vertex &vs0 = s.new_element<Vertex>(2.0, 1.0);
        Vertex &ve0 = s.new_element<Vertex>(1.0, 2.0);

        CircularArc &c0 = s.new_element<CircularArc>(vs0, ve0, vc0, 1.0);

        EXPECT_TRUE(c.is_identical(&c0, &vc0, 180.0));
        EXPECT_TRUE(c.is_identical(&c0, &vc0, -180.0));
        EXPECT_FALSE(c.is_identical(&c0, &vc0, 179.0));
        EXPECT_FALSE(c.is_identical(&c0, &vc0, 181.0));

        Vertex &v1origin = s.new_element<Vertex>(2.0, 2.0);
        Vertex &vc1 = s.new_element<Vertex>(1.0, 3.0);
        Vertex &vs1 = s.new_element<Vertex>(1.0, 4.0);
        Vertex &ve1 = s.new_element<Vertex>(0.0, 3.0);

        CircularArc &c1 = s.new_element<CircularArc>(vs1, ve1, vc1, 1.0);

        EXPECT_TRUE(c.is_identical(&c1, &v1origin, 90.0));
        EXPECT_TRUE(c.is_identical(&c1, &v1origin, -270.0));
        EXPECT_FALSE(c.is_identical(&c1, &v1origin, 89.0));
        EXPECT_FALSE(c.is_identical(&c1, &v1origin, 91.0));

        CircularArc &c2 = s.new_element<CircularArc>(ve0, vs0, vc0, 1.0); // reverse

        EXPECT_FALSE(c.is_identical(&c2, &vc0, 180.0));
        EXPECT_FALSE(c.is_identical(&c2, &vc0, -180.0));

        CircularArc &c3 = s.new_element<CircularArc>(ve1, vs1, vc1, 1.0);

        EXPECT_FALSE(c.is_identical(&c3, &v1origin, 90.0));
        EXPECT_FALSE(c.is_identical(&c3, &v1origin, -270.0));
    }
}

TEST(CIRCULARARC, METHOD_is_coincident) {
    {   // ARGS::(CircularArc)
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

        CircularArc &c0 = s.new_element<CircularArc>(v5, v1, v0, M_SQRT2);
        CircularArc &c1 = s.new_element<CircularArc>(v6, v1, v0, M_SQRT2);
        CircularArc &c2 = s.new_element<CircularArc>(v7, v3, v0, 2.0 * M_SQRT2);
        CircularArc &c3 = s.new_element<CircularArc>(v5, v0, v8, M_SQRT2);

        EXPECT_TRUE(c0.is_coincident(&c1));
        EXPECT_FALSE(c0.is_coincident(&c2));
        EXPECT_FALSE(c0.is_coincident(&c3));
    }

    {   // ARGS::(LineSegment)
        Sketch s;

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(1.0, 1.0);
        Vertex &v5 = s.new_element<Vertex>(M_SQRT2, 0.0);

        CircularArc &c0 = s.new_element<CircularArc>(v5, v1, v0, M_SQRT2);
        LineSegment l0 = LineSegment();

        EXPECT_FALSE(c0.is_coincident(&l0));
    }
}