#include <math.h>
#include <tgmath.h>
#include "test_Sketch.hpp"

TEST(CircularArc, constructor) {
    { //ARGS::()
        EXPECT_NO_THROW(CircularArc c);
    }

    { //ARGS::(Vertex,Vertex,Vertex)
        std::shared_ptr<Vertex> v0, v1, center;
        EXPECT_NO_THROW(CircularArc c(v0, v1, center));
    }
}

TEST(CircularArc, point) {
    {
        double r = 1.0;
        double a0 = -M_PI_2;
        double a1 = M_PI_2;

        std::shared_ptr<Vertex> vc = std::make_shared<Vertex>(0.0, 0.0);
        std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(r * cos(a0), r * sin(a0));
        std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(r * cos(a1), r * sin(a1));
        CircularArc c{v0, v1, vc, r};

        for (size_t i = 0; i < 11; ++i) {
            double s = i / 10.0;
            double a = a0 * (1.0 - s) + a1 * s;

            sPoint v = c.point(s);
            EXPECT_NEAR(r * cos(a), v.x(), TOL);
            EXPECT_NEAR(r * sin(a), v.y(), TOL);
        }
    }

    {
        double r = 1.0;
        double a0 = M_PI_2;
        double a1 = M_PI_2 + M_PI;

        std::shared_ptr<Vertex> vc = std::make_shared<Vertex>(0.0, 0.0);
        std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(r * cos(a0), r * sin(a0));
        std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(r * cos(a1), r * sin(a1));
        CircularArc c{v0, v1, vc, r};

        for (size_t i = 0; i < 11; ++i) {
            double s = i / 10.0;
            double a = a0 * (1.0 - s) + a1 * s;

            sPoint v = c.point(s);
            EXPECT_NEAR(r * cos(a), v.x(), TOL);
            EXPECT_NEAR(r * sin(a), v.y(), TOL);
        }
    }
}

TEST(CircularArc, tangent) {
    {
        double r = 1.0;
        double a0 = -M_PI_2;
        double a1 = M_PI_2;

        std::shared_ptr<Vertex> vc = std::make_shared<Vertex>(0.0, 0.0);
        std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(r * cos(a0), r * sin(a0));
        std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(r * cos(a1), r * sin(a1));
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

        std::shared_ptr<Vertex> vc = std::make_shared<Vertex>(0.0, 0.0);
        std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(r * cos(a0), r * sin(a0));
        std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(r * cos(a1), r * sin(a1));
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

TEST(CircularArc, a) {
    {
        double r = 1.0;
        double a0 = -M_PI_2;
        double a1 = M_PI_2;

        std::shared_ptr<Vertex> vc = std::make_shared<Vertex>(0.0, 0.0);
        std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(r * cos(a0), r * sin(a0));
        std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(r * cos(a1), r * sin(a1));
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

        std::shared_ptr<Vertex> vc = std::make_shared<Vertex>(0.0, 0.0);
        std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(r * cos(a0), r * sin(a0));
        std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(r * cos(a1), r * sin(a1));
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

TEST(CircularArc, da) {
    {
        double r = 1.0;
        double a0 = -M_PI_2;
        double a1 = M_PI_2;

        std::shared_ptr<Vertex> vc = std::make_shared<Vertex>(0.0, 0.0);
        std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(r * cos(a0), r * sin(a0));
        std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(r * cos(a1), r * sin(a1));
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

        std::shared_ptr<Vertex> vc = std::make_shared<Vertex>(0.0, 0.0);
        std::shared_ptr<Vertex> v0 = std::make_shared<Vertex>(r * cos(a0), r * sin(a0));
        std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(r * cos(a1), r * sin(a1));
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

TEST(CircularArc, on_manifold) {
    { //ARGS::(Vertex)
        Sketch s;

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(1.0, 1.0);
        auto v2 = s.new_element<Vertex>(1.1, 1.1);
        auto v3 = s.new_element<Vertex>(M_SQRT2, 0.0);
        auto v4 = s.new_element<Vertex>(0.0, M_SQRT2);
        auto v5 = s.new_element<Vertex>(0.5, 0.5);
        auto v6 = s.new_element<Vertex>(0.5, sqrt(2.0 - 0.25));

        auto c0 = s.new_element<CircularArc>(v3, v1, v0, M_SQRT2);

        EXPECT_FALSE(c0->on_manifold(v0));
        EXPECT_TRUE(c0->on_manifold(v1));
        EXPECT_FALSE(c0->on_manifold(v2));
        EXPECT_TRUE(c0->on_manifold(v3));
        EXPECT_TRUE(c0->on_manifold(v4));
        EXPECT_FALSE(c0->on_manifold(v5));
        EXPECT_TRUE(c0->on_manifold(v6));

        s.delete_me();
    }

    { //ARGS::(Vertex,Vertex,double)
        Sketch s;

        auto vcc = s.new_element<Vertex>(0.0, 0.0);
        auto vc0 = s.new_element<Vertex>(1.0, 0.0);
        auto vc1 = s.new_element<Vertex>(0.0, 1.0);
        double rc = 1.0;

        auto c = s.new_element<CircularArc>(vc0, vc1, vcc, rc);

        auto v0 = s.new_element<Vertex>(M_SQRT1_2, M_SQRT1_2);
        double a0 = 90.0;
        EXPECT_TRUE(c->on_manifold(v0, vcc, a0));
        EXPECT_TRUE(c->on_manifold(v0, vcc, a0 + 1.0));
        EXPECT_TRUE(c->on_manifold(v0, vcc, a0 - 1.0));

        double r1 = hypot(M_SQRT1_2 - 1.0, M_SQRT1_2);
        double a1 = atan2(M_SQRT1_2, M_SQRT1_2 - 1.0) * 180.0 / M_PI;
        auto v1 = s.new_element<Vertex>(1.0 + r1, 0.0);

        EXPECT_TRUE(c->on_manifold(v1, vc0, a1));
        EXPECT_FALSE(c->on_manifold(v1, vc0, a1 + 1.0));
        EXPECT_FALSE(c->on_manifold(v1, vc0, a1 - 1.0));

        auto v2 = s.new_element<Vertex>(0.0, 3.0);
        auto vo2 = s.new_element<Vertex>(0.0, 2.0);
        double a2 = -180.0;

        EXPECT_TRUE(c->on_manifold(v2, vo2, a2));
        EXPECT_FALSE(c->on_manifold(v2, vo2, a2 + 1.0));
        EXPECT_FALSE(c->on_manifold(v2, vo2, a2 - 1.0));

        auto v3 = s.new_element<Vertex>(3.0, 0.0);
        auto vo3 = s.new_element<Vertex>(2.0, 0.0);
        double a3 = 180.0;

        EXPECT_TRUE(c->on_manifold(v3, vo3, a3));
        EXPECT_FALSE(c->on_manifold(v3, vo3, a3 + 1.0));
        EXPECT_FALSE(c->on_manifold(v3, vo3, a3 - 1.0));

        s.delete_me();
    }
}

TEST(CircularArc, on_segment) {
    { //ARGS::(Vertex)
        Sketch s;

        double r = 1.0;
        auto v0 = s.new_element<Vertex>(1.0, 2.0);
        auto v1 = s.new_element<Vertex>(1.0 + r, 2.0);
        auto v2 = s.new_element<Vertex>(1.0, 2.0 + r);

        auto c = s.new_element<CircularArc>(v1, v2, v0, r);

        auto von1 = v1;
        auto von2 = v2;
        auto von3 = s.new_element<Vertex>(1.0 + r * cos(M_PI_4), 2.0 + r * sin(M_PI_4));

        auto voff1 = v0;
        auto voff2 = s.new_element<Vertex>(1.0 + r * cos(-M_PI_4), 2.0 + r * sin(-M_PI_4));
        auto voff3 = s.new_element<Vertex>(1.0 + r * cos(3.0 * M_PI_4), 2.0 + r * sin(3.0 * M_PI_4));
        auto voff4 = s.new_element<Vertex>(1.0 + 2.0 * r * cos(M_PI_4), 2.0 + 2.0 * r * sin(M_PI_4));

        EXPECT_TRUE(c->on_segment(von1));
        EXPECT_TRUE(c->on_segment(von2));
        EXPECT_TRUE(c->on_segment(von3));

        EXPECT_FALSE(c->on_segment(voff1));
        EXPECT_FALSE(c->on_segment(voff2));
        EXPECT_FALSE(c->on_segment(voff3));
        EXPECT_FALSE(c->on_segment(voff4));

        s.delete_me();
    }

    {   // ARGS::(Vertex,Vertex,double)
        Sketch s;

        double r = 1.0;
        auto v0 = s.new_element<Vertex>(1.0, 2.0);
        auto v1 = s.new_element<Vertex>(1.0 + r, 2.0);
        auto v2 = s.new_element<Vertex>(1.0, 2.0 + r);

        auto c = s.new_element<CircularArc>(v1, v2, v0, r);

        auto origin = s.new_element<Vertex>(0.0, 0.0);

        auto von1 = s.new_element<Vertex>(1.0 + r, -2.0);
        double aon1 = 2.0 * atan2(-von1->y(), von1->x()) * 180.0 / M_PI;

        auto von2 = s.new_element<Vertex>(1.0, -2.0 - r);
        double aon2 = 2.0 * atan2(-von2->y(), von2->x()) * 180.0 / M_PI;

        auto von3 = s.new_element<Vertex>(1.0 + r * cos(M_PI_4), -2.0 - r * sin(M_PI_4));
        double aon3 = 2.0 * atan2(-von3->y(), von3->x()) * 180.0 / M_PI;

        auto voff1 = s.new_element<Vertex>(1.0, -2.0);
        double aoff1 = 2.0 * atan2(-voff1->y(), voff1->x()) * 180.0 / M_PI;

        auto voff2 = s.new_element<Vertex>(1.0 + r * cos(-M_PI_4), -2.0 - r * sin(-M_PI_4));
        double aoff2 = 2.0 * atan2(-voff2->y(), voff2->x()) * 180.0 / M_PI;

        auto voff3 = s.new_element<Vertex>(1.0 + r * cos(3.0 * M_PI_4), -2.0 - r * sin(3.0 * M_PI_4));
        double aoff3 = 2.0 * atan2(-voff3->y(), voff3->x()) * 180.0 / M_PI;

        auto voff4 = s.new_element<Vertex>(1.0 + 2.0 * r * cos(M_PI_4), -2.0 - 2.0 * r * sin(M_PI_4));
        double aoff4 = 2.0 * atan2(-voff4->y(), voff4->x()) * 180.0 / M_PI;

        EXPECT_TRUE(c->on_segment(von1, origin, aon1));
        EXPECT_TRUE(c->on_segment(von2, origin, aon2));
        EXPECT_TRUE(c->on_segment(von3, origin, aon3));

        EXPECT_FALSE(c->on_segment(voff1, origin, aoff1));
        EXPECT_FALSE(c->on_segment(voff2, origin, aoff2));
        EXPECT_FALSE(c->on_segment(voff3, origin, aoff3));
        EXPECT_FALSE(c->on_segment(voff4, origin, aoff4));

        s.delete_me();
    }
}

TEST(CircularArc, is_identical) {
    {   //ARGS::(Vertex)
        Sketch s;

        auto vcc = s.new_element<Vertex>(0.0, 0.0);
        auto vc0 = s.new_element<Vertex>(1.0, 0.0);
        auto vc1 = s.new_element<Vertex>(0.0, 1.0);

        auto c = s.new_element<CircularArc>(vc0, vc1, vcc, 1.0);

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(1.0, 0.0);
        auto v2 = s.new_element<Vertex>(0.0, 1.0);
        auto v3 = s.new_element<Vertex>(1.0, 1.0);
        auto v4 = s.new_element<Vertex>(M_SQRT1_2, M_SQRT1_2);

        auto c0 = s.new_element<CircularArc>(v1, v2, v0, 1.0);

        // True
        EXPECT_TRUE(c->is_identical(c));
        EXPECT_TRUE(c->is_identical(c0));

        // #TODO:	Radius does not match distance of endpoints from center
        //			Could be identical depending on other constraints
        //			Behavior is undefined unless Sketch::solve() is called
        auto c1 = s.new_element<CircularArc>(vc0, vc1, vcc, 0.5);
        EXPECT_FALSE(c->is_identical(c1));

        auto c2 = s.new_element<CircularArc>(v1, v2, v0, 0.5);
        EXPECT_FALSE(c->is_identical(c2));

        // False
        auto c3 = s.new_element<CircularArc>(vc1, vc0, vcc, 0.5);
        auto c4 = s.new_element<CircularArc>(v2, v1, v0, 1.0);
        auto c5 = s.new_element<CircularArc>(v1, v2, v3, 1.0);
        auto c6 = s.new_element<CircularArc>(v2, v1, v3, 1.0);

        EXPECT_FALSE(c->is_identical(c3));
        EXPECT_FALSE(c->is_identical(c4));
        EXPECT_FALSE(c->is_identical(c5));
        EXPECT_FALSE(c->is_identical(c6));

        s.delete_me();
    }

    {   // ARGS::(Vertex,Vertex,double)
        Sketch s;

        auto vc = s.new_element<Vertex>(1.0, 1.0);
        auto vs = s.new_element<Vertex>(0.0, 1.0);
        auto ve = s.new_element<Vertex>(1.0, 0.0);

        auto c = s.new_element<CircularArc>(vs, ve, vc, 1.0);

        auto vc0 = s.new_element<Vertex>(1.0, 1.0);
        auto vs0 = s.new_element<Vertex>(2.0, 1.0);
        auto ve0 = s.new_element<Vertex>(1.0, 2.0);

        auto c0 = s.new_element<CircularArc>(vs0, ve0, vc0, 1.0);

        EXPECT_TRUE(c->is_identical(c0, vc0, 180.0));
        EXPECT_TRUE(c->is_identical(c0, vc0, -180.0));
        EXPECT_FALSE(c->is_identical(c0, vc0, 179.0));
        EXPECT_FALSE(c->is_identical(c0, vc0, 181.0));

        auto v1origin = s.new_element<Vertex>(2.0, 2.0);
        auto vc1 = s.new_element<Vertex>(1.0, 3.0);
        auto vs1 = s.new_element<Vertex>(1.0, 4.0);
        auto ve1 = s.new_element<Vertex>(0.0, 3.0);

        auto c1 = s.new_element<CircularArc>(vs1, ve1, vc1, 1.0);

        EXPECT_TRUE(c->is_identical(c1, v1origin, 90.0));
        EXPECT_TRUE(c->is_identical(c1, v1origin, -270.0));
        EXPECT_FALSE(c->is_identical(c1, v1origin, 89.0));
        EXPECT_FALSE(c->is_identical(c1, v1origin, 91.0));

        auto c2 = s.new_element<CircularArc>(ve0, vs0, vc0, 1.0); // reverse

        EXPECT_FALSE(c->is_identical(c2, vc0, 180.0));
        EXPECT_FALSE(c->is_identical(c2, vc0, -180.0));

        auto c3 = s.new_element<CircularArc>(ve1, vs1, vc1, 1.0);

        EXPECT_FALSE(c->is_identical(c3, v1origin, 90.0));
        EXPECT_FALSE(c->is_identical(c3, v1origin, -270.0));

        s.delete_me();
    }
}

TEST(CircularArc, is_coincident) {
    {   // ARGS::(CircularArc)
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

        auto c0 = s.new_element<CircularArc>(v5, v1, v0, M_SQRT2);
        auto c1 = s.new_element<CircularArc>(v6, v1, v0, M_SQRT2);
        auto c2 = s.new_element<CircularArc>(v7, v3, v0, 2.0 * M_SQRT2);
        auto c3 = s.new_element<CircularArc>(v5, v0, v8, M_SQRT2);

        EXPECT_TRUE(c0->is_coincident(c1));
        EXPECT_FALSE(c0->is_coincident(c2));
        EXPECT_FALSE(c0->is_coincident(c3));

        s.delete_me();
    }

    {   // ARGS::(LineSegment)
        Sketch s;

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(1.0, 1.0);
        auto v5 = s.new_element<Vertex>(M_SQRT2, 0.0);

        auto c0 = s.new_element<CircularArc>(v5, v1, v0, M_SQRT2);
        auto l0 = std::make_shared<LineSegment>();

        EXPECT_FALSE(c0->is_coincident(l0));

        s.delete_me();
    }
}