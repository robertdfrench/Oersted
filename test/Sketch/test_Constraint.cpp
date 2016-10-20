#include "test_Sketch.hpp"

TEST(Constraint, Vertex) {
    Sketch s;

    Vertex &v = s.new_element<Vertex>(M_E, M_PI);

    Fixation &f = s.new_element<Fixation>(v);

    s.solve();

    EXPECT_NEAR(M_E, v.x(), TOL);
    EXPECT_NEAR(M_PI, v.y(), TOL);
}

TEST(Constraint, Fixation_Vertex_Length_Line) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(M_1_PI, M_2_SQRTPI);
    Vertex &v1 = s.new_element<Vertex>(M_E, M_LN10);

    LineSegment &line = s.new_element<LineSegment>(v0, v1);

    Fixation &fix = s.new_element<Fixation>(v0);
    Length &len = s.new_element<Length>(line, M_LN2);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Fixation_Vertex_Length_Line.csv");

    EXPECT_NEAR(M_1_PI, v0.x(), TOL);
    EXPECT_NEAR(M_2_SQRTPI, v0.y(), TOL);
    EXPECT_NEAR(M_LN2, line.length(), TOL);
}

TEST(Constraint, Line_Horizontal) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(M_1_PI, M_2_SQRTPI);
    Vertex &v1 = s.new_element<Vertex>(M_E, M_LN10);

    LineSegment &line = s.new_element<LineSegment>(v0, v1);
    Horizontal &horz = s.new_element<Horizontal>(line);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Line_Horizontal.csv");

    EXPECT_NEAR(v0.y(), v1.y(), TOL);
}

TEST(Constraint, Line_Vertical) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(M_1_PI, M_2_SQRTPI);
    Vertex &v1 = s.new_element<Vertex>(M_E, M_LN10);

    LineSegment &line = s.new_element<LineSegment>(v0, v1);

    Vertical &vert = s.new_element<Vertical>(line);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Vertical.csv");

    EXPECT_NEAR(v0.x(), v1.x(), TOL);
}

TEST(Constraint, Line_Length) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(3.14159, 2.7183);
    Vertex &v1 = s.new_element<Vertex>(6.14159, 6.7183);

    LineSegment &line = s.new_element<LineSegment>(v0, v1);

    Length &length = s.new_element<Length>(line, 1.0);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Line_Length.csv");

    EXPECT_NEAR(1.0, line.length(), TOL);
}

TEST(Constraint, internal_CircularArc) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.1, 0.9);
    Vertex &v1 = s.new_element<Vertex>(0.8, 0.2);
    Vertex &vc = s.new_element<Vertex>(0.4, 0.5);

    CircularArc &c = s.new_element<CircularArc>(v0, v1, vc, 2.1);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__internal_CircularArc.csv");

    EXPECT_NEAR(c.radius(), hypot(v0.x() - vc.x(), v0.y() - vc.y()), TOL);
    EXPECT_NEAR(c.radius(), hypot(v1.x() - vc.x(), v1.y() - vc.y()), TOL);
}

TEST(Constraint, Radius_CircularArc) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.1, 0.9);
    Vertex &v1 = s.new_element<Vertex>(0.8, 0.2);
    Vertex &vc = s.new_element<Vertex>(0.4, 0.5);

    CircularArc &c = s.new_element<CircularArc>(v0, v1, vc, 2.1);

    Radius &r = s.new_element<Radius>(c, 3.14);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Radius_CircularArc.csv");

    EXPECT_NEAR(c.radius(), hypot(v0.x() - vc.x(), v0.y() - vc.y()), TOL);
    EXPECT_NEAR(c.radius(), hypot(v1.x() - vc.x(), v1.y() - vc.y()), TOL);
    EXPECT_NEAR(3.14, c.radius(), TOL);
    EXPECT_NEAR(3.14, hypot(v0.x() - vc.x(), v0.y() - vc.y()), TOL);
    EXPECT_NEAR(3.14, hypot(v1.x() - vc.x(), v1.y() - vc.y()), TOL);
}

TEST(Constraint, Tangency_CircularArc_LineSegment_0) {
    // Test when CircularArc and LineSegment have no shared endpoints
    Sketch s;

    Vertex &vc0 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &vc1 = s.new_element<Vertex>(0.0, 1.0);
    Vertex &vcc = s.new_element<Vertex>(0.0, 0.0);

    CircularArc &c = s.new_element<CircularArc>(vc0, vc1, vcc, 1.0);

    Fixation &fix = s.new_element<Fixation>(vcc);
    Radius &rad = s.new_element<Radius>(c, 1.0);

    Vertex &vv0 = s.new_element<Vertex>(M_1_PI, M_2_SQRTPI);
    Vertex &vv1 = s.new_element<Vertex>(M_E, M_LN10);

    LineSegment &lv = s.new_element<LineSegment>(vv0, vv1);

    Vertical &vert = s.new_element<Vertical>(lv);
    Tangency &tanv = s.new_element<Tangency>(c, lv);

    Vertex &vh0 = s.new_element<Vertex>(M_1_PI, M_2_SQRTPI);
    Vertex &vh1 = s.new_element<Vertex>(M_E, M_LN10);

    LineSegment &lh = s.new_element<LineSegment>(vh0, vh1);

    Horizontal &horz = s.new_element<Horizontal>(lh);
    Tangency &tanh = s.new_element<Tangency>(c, lh);

    Vertex &v450 = s.new_element<Vertex>(M_SQRT2, 0.0);
    Vertex &v451 = s.new_element<Vertex>(-M_1_PI, M_SQRT2);

    Fixation &fix45 = s.new_element<Fixation>(v450);

    LineSegment &l45 = s.new_element<LineSegment>(v450, v451);
    LineSegment &l45v = s.new_element<LineSegment>(vcc, v451);

    Tangency &tan45 = s.new_element<Tangency>(c, l45);
    Vertical &vert45 = s.new_element<Vertical>(l45v);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Tangency_CircularArc_LineSegment_0.csv");

    EXPECT_NEAR(0.0, vcc.x(), TOL);
    EXPECT_NEAR(0.0, vcc.y(), TOL);
    EXPECT_NEAR(1.0, c.radius(), TOL);

    EXPECT_NEAR(vv0.x(), vv1.x(), TOL);
    EXPECT_NEAR(vh0.y(), vh1.y(), TOL);
    EXPECT_NEAR(M_SQRT2, v450.x(), TOL);
    EXPECT_NEAR(0.0, v450.y(), TOL);

    EXPECT_NEAR(1.0, vv0.x(), TOL);
    EXPECT_NEAR(1.0, vv1.x(), TOL);

    EXPECT_NEAR(1.0, vh0.y(), TOL);
    EXPECT_NEAR(1.0, vh1.y(), TOL);

    EXPECT_NEAR((M_SQRT2 - v451.x()), v451.y(), TOL);
    EXPECT_NEAR(0.0, v451.x(), TOL);
    EXPECT_NEAR(M_SQRT2, v451.y(), TOL);
}

TEST(Constraint, Tangency_CircularArc_LineSegment_1) {
    // Test when CircularArc and LineSegment have shared endpoints
    Sketch s;

    Vertex &vc = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v0 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(0.0, 1.0);
    Vertex &v2 = s.new_element<Vertex>(M_SQRT2, 0.0);

    CircularArc &c = s.new_element<CircularArc>(v0, v1, vc, 1.0);
    LineSegment &l = s.new_element<LineSegment>(v1, v2);

    Fixation &fc = s.new_element<Fixation>(vc);
    Fixation &f0 = s.new_element<Fixation>(v2);
    Tangency &t = s.new_element<Tangency>(c, l);
    Radius &r = s.new_element<Radius>(c, 1.0);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__CircularArc_LineSegment_1.csv");

    EXPECT_NEAR(1.0, c.radius(), TOL);
    EXPECT_NEAR(M_PI_2, (atan2(v1.y() - vc.y(), v1.x() - vc.x()) - atan2(v2.y() - v1.y(), v2.x() - v1.x())), TOL * M_PI_2);
    EXPECT_NEAR(M_SQRT1_2, v1.x(), TOL * M_SQRT1_2);
    EXPECT_NEAR(M_SQRT1_2, v1.y(), TOL * M_SQRT1_2);
}

TEST(Constraint, Angle_LineSegment_LineSegment) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v00 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v01 = s.new_element<Vertex>(0.0, 1.0);

    LineSegment &line0 = s.new_element<LineSegment>(v0, v00);
    LineSegment &line1 = s.new_element<LineSegment>(v0, v01);

    Fixation &fix = s.new_element<Fixation>(v0);
    Horizontal &horz = s.new_element<Horizontal>(line0);
    Angle &ang = s.new_element<Angle>(line0, line1, 30.0);
    Length &len = s.new_element<Length>(line1, 1.0);

    for (double i = 0.0; i < 16.0; i++) {
        double a = 22.5 * i;

        ang.Dim = a;

        s.solve();

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint_Angle_LineSegment_LineSegment.csv");

        EXPECT_NEAR(0.0, v0.x(), TOL);
        EXPECT_NEAR(0.0, v0.y(), TOL);
        EXPECT_NEAR(v0.y(), v00.y(), TOL);
        EXPECT_NEAR(1.0, line1.length(), TOL);

        EXPECT_NEAR(tan(M_PI * i / 8.0), ((v01.y() - v0.y()) / (v01.x() - v0.x())), TOL);
        EXPECT_NEAR(cos(M_PI * i / 8.0), v01.x(), TOL);
        EXPECT_NEAR(sin(M_PI * i / 8.0), v01.y(), TOL);
    }
}

TEST(Constraint, Angle_Coincident_LineSegment) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);

    Vertex &v00 = s.new_element<Vertex>(0.5, 0.0);
    Vertex &v01 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v10 = s.new_element<Vertex>(0.5, 0.5);
    Vertex &v11 = s.new_element<Vertex>(1.0, 1.0);

    LineSegment &line0 = s.new_element<LineSegment>(v00, v01);
    LineSegment &line1 = s.new_element<LineSegment>(v10, v11);

    Fixation &fix = s.new_element<Fixation>(v0);
    Horizontal &horz = s.new_element<Horizontal>(line0);
    Angle &ang = s.new_element<Angle>(line0, line1, 0.0);
    Length &len = s.new_element<Length>(line1, 0.5);
    Coincident<LineSegment> &coin0 = s.new_element<Coincident<LineSegment>>(v0, line0);
    Coincident<LineSegment> &coin1 = s.new_element<Coincident<LineSegment>>(v0, line1);

    size_t N = 15;
    double dar = 2.0 * M_PI / N;
    double dad = 360.0 / N;
    for (size_t i = 0; i < N; ++i) {
        double a = i * dad;

        ang.Dim = a;

        s.solve();

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Angle_Coincident.csv");

        EXPECT_NEAR(0.0, v0.x(), TOL);
        EXPECT_NEAR(0.0, v0.y(), TOL);
        EXPECT_NEAR(v00.y(), v01.y(), TOL);
        EXPECT_NEAR(0.5, line1.length(), TOL);

        EXPECT_NEAR(tan(dar * i) * (v11.x() - v10.x()), (v11.y() - v10.y()), TOL);
        EXPECT_NEAR(0.0, tan(dar * i) * v11.x() - v11.y(), TOL);
        EXPECT_NEAR(0.0, tan(dar * i) * v10.x() - v10.y(), TOL);
    }
}

TEST(Constraint, Distance_Vertex_Vertex) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(M_1_PI, M_2_SQRTPI);
    Vertex &v1 = s.new_element<Vertex>(M_E, M_LN10);

    Distance<Vertex> &d = s.new_element<Distance<Vertex>>(v0, v1, M_LN2);

    s.solve();

    EXPECT_NEAR(M_LN2, hypot(v0.x() - v1.x(), v0.y() - v1.y()), TOL * M_LN2);
}

TEST(Constraint, Distance_LineSegment) {
    Sketch s;

    Vertex &v00 = s.new_element<Vertex>(M_1_PI, M_LOG10E);
    Vertex &v01 = s.new_element<Vertex>(M_2_SQRTPI, M_LOG10E);
    Vertex &v10 = s.new_element<Vertex>(M_E, M_LOG2E);
    Vertex &v11 = s.new_element<Vertex>(M_LN2, M_LOG2E);

    LineSegment &l0 = s.new_element<LineSegment>(v00, v01);
    LineSegment &l1 = s.new_element<LineSegment>(v10, v11);
    Distance<LineSegment> &d = s.new_element<Distance<LineSegment>>(l0, l1, M_PI);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Distance_LineSegment.csv");

    double v0x = v01.x() - v00.x();
    double v1x = v11.x() - v10.x();
    double v0y = v01.y() - v00.y();
    double v1y = v11.y() - v10.y();

    double d0 = sqrt(v0x * v0x + v0y * v0y);
    v0x /= d0;
    v0y /= d0;

    double d1 = sqrt(v1x * v1x + v1y * v1y);
    v1x /= d1;
    v1y /= d1;

    double dot = abs(v0x * v1x + v0y * v1y);
    double cross = abs(v0x * v1y - v0y * v1x);

    double len;
    len = ((v01.x() - v00.x()) * (v11.y() - v00.y()) - (v01.y() - v00.y()) * (v11.x() - v00.x()));
    len -= ((v11.x() - v00.x()) * (v10.y() - v00.y()) - (v11.y() - v00.y()) * (v10.x() - v00.x()));
    len = len / (d0 + d1);

    EXPECT_NEAR(1.0, dot, TOL);
    EXPECT_NEAR(0.0, cross, TOL);
    EXPECT_NEAR(M_PI, len, TOL * M_PI);
}

TEST(Constraint, Distance_LineSegment_initially_perpendicular) {
    Sketch s;

    Vertex &v00 = s.new_element<Vertex>(-1.0, 0.0);
    Vertex &v01 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v10 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v11 = s.new_element<Vertex>(1.0, 1.0);

    LineSegment &l0 = s.new_element<LineSegment>(v00, v01);
    LineSegment &l1 = s.new_element<LineSegment>(v10, v11);
    Distance<LineSegment> &d = s.new_element<Distance<LineSegment>>(l0, l1, 1.0);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Distance_LineSegment_initially_perpendicular");

    double dx0 = v01.x() - v00.x();
    double dy0 = v01.y() - v00.y();
    double dr0 = sqrt(dx0 * dx0 + dy0 * dy0);
    dx0 /= dr0;
    dy0 /= dr0;

    double dx1 = v11.x() - v10.x();
    double dy1 = v11.y() - v10.y();
    double dr1 = sqrt(dx1 * dx1 + dy1 * dy1);
    dx1 /= dr1;
    dy1 /= dr1;

    double dot = abs(dx0 * dx1 + dy0 * dy1);

    EXPECT_NEAR(1.0, dot, TOL);

    double x, y, l;
    x = (v01.x() + v00.x()) / 2.0;
    y = (v01.y() + v00.y()) / 2.0;

    dx0 = (v10.x() - x);
    dx1 = (v11.x() - x);
    dy0 = (v10.y() - y);
    dy1 = (v11.y() - y);

    l = abs(dx0 * dy1 - dy0 * dx1) / dr1;
    EXPECT_NEAR(1.0, l, TOL);

    x = (v11.x() + v10.x()) / 2.0;
    y = (v11.y() + v10.y()) / 2.0;

    dx0 = (v00.x() - x);
    dx1 = (v01.x() - x);
    dy0 = (v00.y() - y);
    dy1 = (v01.y() - y);

    l = abs(dx0 * dy1 - dy0 * dx1) / dr0;
    EXPECT_NEAR(1.0, l, TOL);
}

TEST(Constraint, Distance_LineSegment_midpoint_intersection) {
    Sketch s;

    Vertex &v00 = s.new_element<Vertex>(-1.0, 0.0);
    Vertex &v01 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v10 = s.new_element<Vertex>(0.0, -1.0);
    Vertex &v11 = s.new_element<Vertex>(0.0, 1.0);

    LineSegment &l0 = s.new_element<LineSegment>(v00, v01);
    LineSegment &l1 = s.new_element<LineSegment>(v10, v11);
    Distance<LineSegment> &d = s.new_element<Distance<LineSegment>>(l0, l1, 1.0);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Distance_LineSegment_midpoint_intersection.csv");

    double dx0 = v01.x() - v00.x();
    double dy0 = v01.y() - v00.y();
    double dr0 = sqrt(dx0 * dx0 + dy0 * dy0);
    dx0 /= dr0;
    dy0 /= dr0;

    double dx1 = v11.x() - v10.x();
    double dy1 = v11.y() - v10.y();
    double dr1 = sqrt(dx1 * dx1 + dy1 * dy1);
    dx1 /= dr1;
    dy1 /= dr1;

    double dot = abs(dx0 * dx1 + dy0 * dy1);

    EXPECT_NEAR(1.0, dot, TOL); // #TODO: Fails because lines intersect at midpoint
}

TEST(Constraint, Distance_CircularArc_exterior) {
    Sketch s;

    Vertex &v00 = s.new_element<Vertex>(0.2, 0.0);
    Vertex &v01 = s.new_element<Vertex>(-0.2, 0.0);
    Vertex &vc0 = s.new_element<Vertex>(0.0, 0.0);

    Vertex &vc1 = s.new_element<Vertex>(1.0, 1.0);
    Vertex &v10 = s.new_element<Vertex>(0.8, 1.0);
    Vertex &v11 = s.new_element<Vertex>(1.2, 1.0);

    CircularArc &c0 = s.new_element<CircularArc>(v00, v01, vc0, 0.2);
    CircularArc &c1 = s.new_element<CircularArc>(v10, v11, vc1, 0.2);

    s.solve();

    Distance<CircularArc> &d = s.new_element<Distance<CircularArc>>(c0, c1, 1.0 / 3.0);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Distance_CircularArc_exterior.csv");

    EXPECT_NEAR(1.0 / 3.0, hypot(vc0.x() - vc1.x(), vc0.y() - vc1.y()) - c0.radius() - c1.radius(), TOL / 3.0);
}

TEST(Constraint, Distance_CircularArc_interior_0) {
    Sketch s;

    Vertex &v00 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v01 = s.new_element<Vertex>(-1.0, 0.0);
    Vertex &vc0 = s.new_element<Vertex>(0.0, 0.0);

    Vertex &vc1 = s.new_element<Vertex>(0.2, 0.2);
    Vertex &v10 = s.new_element<Vertex>(0.4, 0.2);
    Vertex &v11 = s.new_element<Vertex>(0.0, 0.2);

    CircularArc &c0 = s.new_element<CircularArc>(v00, v01, vc0, 1.0);
    CircularArc &c1 = s.new_element<CircularArc>(v10, v11, vc1, 0.2);

    s.solve();

    Distance<CircularArc> &d = s.new_element<Distance<CircularArc> >(c0, c1, 1.0 / 3.0);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Distance_CircularArc_interior_0.csv");

    EXPECT_NEAR(1.0 / 3.0, c0.radius() - hypot(vc0.x() - vc1.x(), vc0.y() - vc1.y()) - c1.radius(), TOL / 3.0);
}

TEST(Constraint, Distance_CircularArc_interior_1) {
    Sketch s;

    Vertex &v00 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v01 = s.new_element<Vertex>(-1.0, 0.0);
    Vertex &vc0 = s.new_element<Vertex>(0.0, 0.0);

    Vertex &vc1 = s.new_element<Vertex>(0.2, 0.2);
    Vertex &v10 = s.new_element<Vertex>(0.4, 0.2);
    Vertex &v11 = s.new_element<Vertex>(0.0, 0.2);

    CircularArc &c0 = s.new_element<CircularArc>(v00, v01, vc0, 1.0);
    CircularArc &c1 = s.new_element<CircularArc>(v10, v11, vc1, 0.2);

    s.solve();

    Distance<CircularArc> &d = s.new_element<Distance<CircularArc>>(c1, c0, 1.0 / 3.0);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Distance_CircularArc_interior_1.csv");

    EXPECT_NEAR(1.0 / 3.0, c0.radius() - hypot(vc0.x() - vc1.x(), vc0.y() - vc1.y()) - c1.radius(), TOL / 3.0);
}

TEST(Constraint, Coincident_CircularArc) {
    Sketch s;

    Vertex &vc = s.new_element<Vertex>(M_1_PI, M_2_SQRTPI);
    Vertex &v0 = s.new_element<Vertex>(M_1_PI, M_2_SQRTPI + M_LN2);
    Vertex &v1 = s.new_element<Vertex>(M_1_PI, M_2_SQRTPI - M_LN2);

    CircularArc &ca = s.new_element<CircularArc>(v0, v1, vc, M_LN2);

    Vertex &vp = s.new_element<Vertex>(M_E, M_LN10);

    Coincident<CircularArc> &co = s.new_element<Coincident<CircularArc>>(vp, ca);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Coincident_CircularArc.csv");

    double dx = vc.x() - vp.x();
    double dy = vc.y() - vp.y();
    double d = sqrt(dx * dx + dy * dy);
    EXPECT_NEAR(0.0, d - ca.radius(), TOL * ca.radius());
}

TEST(Constraint, Coincident_LineSegment) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(1.0, 1.0);
    Vertex &vp = s.new_element<Vertex>(1.1, 0.9);

    LineSegment &l = s.new_element<LineSegment>(v0, v1);

    Coincident<LineSegment> &c = s.new_element<Coincident<LineSegment>>(vp, l);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Coincident_LineSegment.csv");

    double dx0 = v0.x() - vp.x();
    double dy0 = v0.y() - vp.y();
    double dx1 = v1.x() - vp.x();
    double dy1 = v1.y() - vp.y();
    double dr0 = sqrt(dx0 * dx0 + dy0 * dy0);
    double dr1 = sqrt(dx1 * dx1 + dy1 * dy1);
    dx0 /= dr0;
    dy0 /= dr0;
    dx1 /= dr1;
    dy1 /= dr1;

    EXPECT_NEAR(0.0, (dx0 * dy1 - dy0 * dx1), TOL);
}

TEST(Constraint, Coincident_LineSegment_deg180) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(1.0, 1.0);
    Vertex &vp = s.new_element<Vertex>(0.25, 0.75);

    LineSegment &l = s.new_element<LineSegment>(v0, v1);

    Coincident<LineSegment> &c = s.new_element<Coincident<LineSegment>>(vp, l);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Constraint__Coincident_LineSegment_deg180.csv");

    double dx0 = v0.x() - vp.x();
    double dy0 = v0.y() - vp.y();
    double dx1 = v1.x() - vp.x();
    double dy1 = v1.y() - vp.y();
    double dr0 = sqrt(dx0 * dx0 + dy0 * dy0);
    double dr1 = sqrt(dx1 * dx1 + dy1 * dy1);
    dx0 /= dr0;
    dy0 /= dr0;
    dx1 /= dr1;
    dy1 /= dr1;

    EXPECT_NEAR(0.0, (dx0 * dy1 - dy0 * dx1), TOL);
}

TEST(Constraint, Parallel) {
    FAIL(); // TODO: Implement test for parallel lines
}

TEST(Constraint, Parallel_initially_perpendicular) {
    FAIL(); // TODO: Implement test for parallel lines that are initially perpendicular
}

TEST(Constraint, Perpendicular) {
    FAIL(); // TODO: Implement test for perpendicular lines
}

TEST(Constraint, Perpendicular_initially_parallel) {
    FAIL(); // TODO: Implement test for perpendicular lines that are initially parallel
}

TEST(Constraint, Symmetry_Horizontal) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v2 = s.new_element<Vertex>(0.41, -0.62);
    Vertex &v3 = s.new_element<Vertex>(0.63, 0.44);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);

    s.new_element<Fixation>(v0);
    s.new_element<Horizontal>(l0);
    s.new_element<Length>(l0, 1.0);

    s.solve();

    EXPECT_NEAR(0.0, v0.x(), TOL);
    EXPECT_NEAR(0.0, v0.y(), TOL);
    EXPECT_NEAR(0.0, v1.y(), TOL);
    EXPECT_NEAR(1.0, l0.length(), TOL);
    EXPECT_NEAR(1.0, v1.x(), TOL);

    s.new_element<Symmetry>(v2, v3, l0);
    s.solve();

    EXPECT_NEAR(v2.x(), v3.x(), TOL);
    EXPECT_NEAR(v2.y(), -v3.y(), TOL);
}

TEST(Constraint, Symmetry_Vertical) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(0.0, 1.0);
    Vertex &v2 = s.new_element<Vertex>(-0.62, 0.41);
    Vertex &v3 = s.new_element<Vertex>(0.44, 0.63);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);

    s.new_element<Fixation>(v0);
    s.new_element<Vertical>(l0);
    s.new_element<Length>(l0, 1.0);

    s.solve();

    EXPECT_NEAR(0.0, v0.x(), TOL);
    EXPECT_NEAR(0.0, v0.y(), TOL);
    EXPECT_NEAR(0.0, v1.x(), TOL);
    EXPECT_NEAR(1.0, l0.length(), TOL);
    EXPECT_NEAR(1.0, v1.y(), TOL);

    s.new_element<Symmetry>(v2, v3, l0);
    s.solve();

    EXPECT_NEAR(v2.x(), -v3.x(), TOL);
    EXPECT_NEAR(v2.y(), v3.y(), TOL);
}

TEST(Constraint, Symmetry_deg45) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(1.0, 1.0);
    Vertex &v2 = s.new_element<Vertex>(0.1, 0.7);
    Vertex &v3 = s.new_element<Vertex>(0.8, 0.2);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    s.new_element<Fixation>(v0);
    s.new_element<Fixation>(v1);

    s.solve();

    EXPECT_NEAR(0.0, v0.x(), TOL);
    EXPECT_NEAR(0.0, v0.y(), TOL);
    EXPECT_NEAR(1.0, v1.x(), TOL);
    EXPECT_NEAR(1.0, v1.y(), TOL);
    EXPECT_NEAR(M_SQRT2, l0.length(), TOL);

    s.new_element<Symmetry>(v2, v3, l0);
    s.solve();

    EXPECT_NEAR(v3.y(), v2.x(), TOL);
    EXPECT_NEAR(v3.x(), v2.y(), TOL);
}

TEST(Constraint, Rotation) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 1.0);
    Vertex &v1 = s.new_element<Vertex>(2.0, 3.0);
    Vertex &v2 = s.new_element<Vertex>(3.0, 4.0);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v0, v2);

    double a = 60.0;
    Rotation &r = s.new_element<Rotation>(v1, v2, v0, a);
    a *= M_PI / 180.0;

    s.solve();

    double dx1 = v1.x() - v0.x();
    double dy1 = v1.y() - v0.y();
    double dx2 = v2.x() - v0.x();
    double dy2 = v2.y() - v0.y();

    double r1 = sqrt(dx1 * dx1 + dy1 * dy1);
    double r2 = sqrt(dx2 * dx2 + dy2 * dy2);

    EXPECT_NEAR(r1, r2, TOL);
    EXPECT_NEAR(v2.x(), dx1 * cos(a) - dy1 * sin(a) + v0.x(), TOL);
    EXPECT_NEAR(v2.y(), dx1 * sin(a) + dy1 * cos(a) + v0.y(), TOL);
}