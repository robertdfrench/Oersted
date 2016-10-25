#include "test_UseCases.hpp"

#define SAVE_DIR "./test/output/UseCases/Rotor/"

TEST(ROTOR, 0) {
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

    Vertex &origin = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v0 = s.new_element<Vertex>(xi, yi);
    Vertex &v1 = s.new_element<Vertex>(xi, -yi);
    Vertex &v2 = s.new_element<Vertex>(xo, -yo);
    Vertex &v3 = s.new_element<Vertex>(xo, yo);

    Vertex &v4 = s.new_element<Vertex>(xo, -ym);
    Vertex &v5 = s.new_element<Vertex>(xo, ym);
    Vertex &v6 = s.new_element<Vertex>(xm, ym);
    Vertex &v7 = s.new_element<Vertex>(xm, -ym);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v3);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    CircularArc &c0 = s.new_element<CircularArc>(v1, v0, origin, ri);
    CircularArc &c1 = s.new_element<CircularArc>(v2, v3, origin, ro);

    s.new_element<LineSegment>(v4, v5);
    s.new_element<LineSegment>(v5, v6);
    s.new_element<LineSegment>(v6, v7);
    s.new_element<LineSegment>(v7, v4);

    s.solve();

    s.build();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "rotor_0d0.csv");

    Mesh m{s};
    m.MaximumElementSize = 2.5;
    m.MinimumElementSize = 0.25;
    m.MinimumElementQuality = M_SQRT1_2;

    m.create();

    m.save_as(SAVE_DIR, "rotor_0_mesh.csv");

    m.refine();

    m.save_as(SAVE_DIR, "rotor_0_mesh_refine.csv");
}

TEST(ROTOR, CIRCULAR_BARRIER_SYNCREL) {
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
    Vertex &origin = s.new_element<Vertex>(0.0, 0.0);

    Vertex &vi0 = s.new_element<Vertex>(ri, 0.0);
    Vertex &vo0 = s.new_element<Vertex>(ro, 0.0);

    Vertex &vi1 = s.new_element<Vertex>(ri * cos(arad), ri * sin(arad));
    Vertex &vo1 = s.new_element<Vertex>(ro * cos(arad), ro * sin(arad));

    Vertex &vm0 = s.new_element<Vertex>(rd0 * cos(arad), rd0 * sin(arad));
    Vertex &vr0 = s.new_element<Vertex>(ro * cos(arad * (1.0 - pd)), ro * sin(arad * (1.0 - pd)));
    Vertex &vc0 = s.new_element<Vertex>(ro * cos(arad), ro * sin(arad));

    Vertex &vm1 = s.new_element<Vertex>(rd1 * cos(arad), rd1 * sin(arad));
    Vertex &vr1 = s.new_element<Vertex>(ro * cos(arad * (1.0 - 2.0 * pd)), ro * sin(arad * (1.0 - 2.0 * pd)));

    // LineSegment
    LineSegment &l0 = s.new_element<LineSegment>(vi0, vo0);

    LineSegment &l1 = s.new_element<LineSegment>(vi1, vm1);
    LineSegment &l2 = s.new_element<LineSegment>(vm1, vm0);
    LineSegment &l3 = s.new_element<LineSegment>(vm0, vo1);

    // CircularArc
    CircularArc &c0 = s.new_element<CircularArc>(vi0, vi1, origin, ri);

    CircularArc &c1 = s.new_element<CircularArc>(vo0, vr1, origin, ro);
    CircularArc &c2 = s.new_element<CircularArc>(vr1, vr0, origin, ro);
    CircularArc &c3 = s.new_element<CircularArc>(vr0, vo1, origin, ro);

    CircularArc &c4 = s.new_element<CircularArc>(vm0, vr0, vc0, pd * dmax);
    CircularArc &c5 = s.new_element<CircularArc>(vm1, vr1, vc0, 2.0 * pd * dmax);

    // Fixation
    Fixation &fo = s.new_element<Fixation>(origin);

    // Horizontal
    Horizontal &ho = s.new_element<Horizontal>(l0);

    // Angles
    Angle &ao = s.new_element<Angle>(l0, l1, adeg);
    Angle &a1 = s.new_element<Angle>(l0, l2, adeg);
    Angle &a2 = s.new_element<Angle>(l0, l3, adeg);

    // Coincident
    Coincident<LineSegment> &coinc0 = s.new_element<Coincident<LineSegment>>(vc0, l1);

    // Radius
    Radius &rad0 = s.new_element<Radius>(c0, ri);
    Radius &rad1 = s.new_element<Radius>(c1, ro);
    Radius &rad2 = s.new_element<Radius>(c2, ro);
    Radius &rad3 = s.new_element<Radius>(c3, ro);

    // Distance
    Distance<Vertex> &dist0 = s.new_element<Distance<Vertex>>(vr0, vo1, pd * dmax);
    Distance<Vertex> &dist1 = s.new_element<Distance<Vertex>>(vr1, vr0, pd * dmax);

    // Solve
    s.solve();

    s.build();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "rotor_1_circular_barrier_syncrel.csv");

    // Mesh
    Mesh m{s};
    m.MaximumElementSize = 1.0;
    m.MinimumElementSize = 0.1;
    m.MinimumElementQuality = M_SQRT1_2;

    m.create();

    m.save_as(SAVE_DIR, "rotor_1_circular_barrier_syncrel_mesh.csv");

    m.refine();

    m.save_as(SAVE_DIR, "rotor_1_circular_barrier_syncrel_mesh_refine.csv");
}