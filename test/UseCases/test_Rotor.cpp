#include "test_UseCases.hpp"

TEST(Rotor, Suite0) {
    Sketch s;

    double ri{25.0};
    double ro{80.0};
    double a{M_PI / 8.0};

    double xi = ri * cos(a);
    double xo = ro * cos(a);
    double yi = ri * sin(a);
    double yo = ro * sin(a);

    double xm = 70.0 * cos(a);
    double ym = yo - (xo - xm);

    auto origin = s.new_element<Vertex>(0.0, 0.0);
    auto v0 = s.new_element<Vertex>(xi, yi);
    auto v1 = s.new_element<Vertex>(xi, -yi);
    auto v2 = s.new_element<Vertex>(xo, -yo);
    auto v3 = s.new_element<Vertex>(xo, yo);

    auto v4 = s.new_element<Vertex>(xo, -ym);
    auto v5 = s.new_element<Vertex>(xo, ym);
    auto v6 = s.new_element<Vertex>(xm, ym);
    auto v7 = s.new_element<Vertex>(xm, -ym);

    auto l0 = s.new_element<LineSegment>(v0, v3);
    auto l1 = s.new_element<LineSegment>(v1, v2);
    auto c0 = s.new_element<CircularArc>(v1, v0, origin, ri);
    auto c1 = s.new_element<CircularArc>(v2, v3, origin, ro);

    s.new_element<LineSegment>(v4, v5);
    s.new_element<LineSegment>(v5, v6);
    s.new_element<LineSegment>(v6, v7);
    s.new_element<LineSegment>(v7, v4);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    s.build();

    s.save_as<SaveMethod::Rasterize>(SDIR, "rotor_0");

    Mesh m{s};
    m.MaximumElementSize = 2.5;
    m.MinimumElementSize = 0.25;
    m.MinimumElementQuality = M_SQRT1_2;

    m.create();

    m.save_as(MDIR, "rotor0");

    m.refine();

    m.save_as(MDIR, "rotor0_refined");
}

TEST(Rotor, Circular_Barrier_Syncrel) {
    Sketch s;

    size_t Np = 8;
    double adeg = 180.0 / Np;
    double arad = adeg * M_PI / 180.0;
    double ri = 25.0;
    double ro = 80.0;
    double dmax = ro * sin(arad);
    double pd = 0.2;
    double rd0 = ro - pd * dmax;
    double rd1 = rd0 - pd * dmax;

    // Vertex
    auto origin = s.new_element<Vertex>(0.0, 0.0);

    auto vi0 = s.new_element<Vertex>(ri, 0.0);
    auto vo0 = s.new_element<Vertex>(ro, 0.0);

    auto vi1 = s.new_element<Vertex>(ri * cos(arad), ri * sin(arad));
    auto vo1 = s.new_element<Vertex>(ro * cos(arad), ro * sin(arad));

    auto vm0 = s.new_element<Vertex>(rd0 * cos(arad), rd0 * sin(arad));
    auto vr0 = s.new_element<Vertex>(ro * cos(arad * (1.0 - pd)), ro * sin(arad * (1.0 - pd)));
    auto vc0 = s.new_element<Vertex>(ro * cos(arad), ro * sin(arad));

    auto vm1 = s.new_element<Vertex>(rd1 * cos(arad), rd1 * sin(arad));
    auto vr1 = s.new_element<Vertex>(ro * cos(arad * (1.0 - 2.0 * pd)), ro * sin(arad * (1.0 - 2.0 * pd)));

    // LineSegment
    auto l0 = s.new_element<LineSegment>(vi0, vo0);

    auto l1 = s.new_element<LineSegment>(vi1, vm1);
    auto l2 = s.new_element<LineSegment>(vm1, vm0);
    auto l3 = s.new_element<LineSegment>(vm0, vo1);

    // CircularArc
    auto c0 = s.new_element<CircularArc>(vi0, vi1, origin, ri);

    auto c1 = s.new_element<CircularArc>(vo0, vr1, origin, ro);
    auto c2 = s.new_element<CircularArc>(vr1, vr0, origin, ro);
    auto c3 = s.new_element<CircularArc>(vr0, vo1, origin, ro);

    auto c4 = s.new_element<CircularArc>(vm0, vr0, vc0, pd * dmax);
    auto c5 = s.new_element<CircularArc>(vm1, vr1, vc0, 2.0 * pd * dmax);

    // Fixation
    auto fo = s.new_element<Fixation>(origin);

    // Horizontal
    auto ho = s.new_element<Horizontal>(l0);

    // Angles
    auto ao = s.new_element<Angle>(l0, l1, adeg);
    auto a1 = s.new_element<Angle>(l0, l2, adeg);
    auto a2 = s.new_element<Angle>(l0, l3, adeg);

    // Coincident
    auto coinc0 = s.new_element<Coincident<LineSegment>>(vc0, l1);

    // Radius
    auto rad0 = s.new_element<Radius>(c0, ri);
    auto rad1 = s.new_element<Radius>(c1, ro);
    auto rad2 = s.new_element<Radius>(c2, ro);
    auto rad3 = s.new_element<Radius>(c3, ro);

    // Distance
    auto dist0 = s.new_element<Distance<Vertex>>(vr0, vo1, pd * dmax);
    auto dist1 = s.new_element<Distance<Vertex>>(vr1, vr0, pd * dmax); // TODO: Fix distance with shared_ptr

    // Solve
    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    s.build();

    s.save_as<SaveMethod::Rasterize>(SDIR, "rotor_circular_barrier_syncrel");

    // Mesh
    Mesh m{s};
    m.MaximumElementSize = 1.0;
    m.MinimumElementSize = 0.1;
    m.MinimumElementQuality = M_SQRT1_2;

    m.create();

    m.save_as(MDIR, "rotor_circular_barrier_syncrel");

    m.refine();

    m.save_as(MDIR, "rotor_circular_barrier_syncrel_refined");
}