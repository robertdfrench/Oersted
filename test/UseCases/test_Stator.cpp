#include "test_UseCases.hpp"

TEST(Stator, 0) {
    Sketch sketch;

    size_t Np{8};
    size_t Nt{3 * Np};

    double ri{80.0};
    double ro{130.0};
    double at{M_PI / Nt};
    double as{at / 2.0};
    double dt{ri * sin(as)};
    double rb{ro - 4.0 * dt};

    auto origin = sketch.new_element_SHARED_PTR<Vertex>(0.0, 0.0);
    auto v0 = sketch.new_element_SHARED_PTR<Vertex>(ri, 0.0);
    auto v1 = sketch.new_element_SHARED_PTR<Vertex>(ro, 0.0);
    auto v2 = sketch.new_element_SHARED_PTR<Vertex>(ro * cos(at), ro * sin(at));
    auto v3 = sketch.new_element_SHARED_PTR<Vertex>(rb * cos(at), rb * sin(at));
    auto v4 = sketch.new_element_SHARED_PTR<Vertex>(ri * cos(at), ri * sin(at));
    auto v5 = sketch.new_element_SHARED_PTR<Vertex>(ri * cos(as), ri * sin(as));
    auto v6 = sketch.new_element_SHARED_PTR<Vertex>(rb * cos(as), rb * sin(as));

    LineSegment &l0 = sketch.new_element<LineSegment>(v0, v1);
    CircularArc &c0 = sketch.new_element<CircularArc>(v1, v2, origin, ro);
    LineSegment &l1 = sketch.new_element<LineSegment>(v2, v3);
    LineSegment &l2 = sketch.new_element<LineSegment>(v3, v4);
    CircularArc &c1 = sketch.new_element<CircularArc>(v5, v4, origin, ri);
    CircularArc &c2 = sketch.new_element<CircularArc>(v0, v5, origin, ri);
    LineSegment &l3 = sketch.new_element<LineSegment>(v5, v6);
    CircularArc &c3 = sketch.new_element<CircularArc>(v6, v3, origin, rb);

    auto mvec = sketch.curves();
    auto m0 = sketch.new_element_SHARED_PTR<MirrorCopy>(mvec, &l1, true);

    auto rvec = sketch.curves();
    auto rcopy = sketch.new_element_SHARED_PTR<RotateCopy>(rvec, origin, 360.0 / Nt, 1, true);

    sketch.solve();

    sketch.save_as<SaveMethod::Rasterize>(SDIR, "stator0_0");

    EXPECT_TRUE(sketch.build());
    EXPECT_EQ(sketch.boundary()->size(), 14);

    Fixation &f = sketch.new_element<Fixation>(origin);

    Coincident<LineSegment> &coin0 = sketch.new_element<Coincident<LineSegment>>(origin, l0);
    Coincident<LineSegment> &coin1 = sketch.new_element<Coincident<LineSegment>>(origin, l1);
    Coincident<LineSegment> &coin2 = sketch.new_element<Coincident<LineSegment>>(origin, l2);

    Horizontal &horz0 = sketch.new_element<Horizontal>(l0);

    Radius &rad0 = sketch.new_element<Radius>(c0, ro);
    Radius &rad1 = sketch.new_element<Radius>(c1, ri);
    Radius &rad2 = sketch.new_element<Radius>(c2, ri);
    Radius &rad3 = sketch.new_element<Radius>(c3, rb);

    Angle &ang0 = sketch.new_element<Angle>(l1, l0, (at * 180.0 / M_PI));
    Angle &ang1 = sketch.new_element<Angle>(l2, l0, (at * 180.0 / M_PI));

    Distance<LineSegment> &dist0 = sketch.new_element<Distance<LineSegment>>(l0, l3, dt);

    sketch.solve();

    sketch.save_as<SaveMethod::Rasterize>(SDIR, "stator0_1");

    EXPECT_TRUE(sketch.build());
    EXPECT_EQ(sketch.boundary()->size(), 14);

    Mesh mesh{sketch};

    mesh.MinimumElementQuality = 0 * M_SQRT1_2;
    mesh.MaximumElementSize = 2.5;
    mesh.MinimumElementSize = 0.25;

    mesh.create();

    EXPECT_TRUE(edges_are_valid(mesh));
    EXPECT_TRUE(edges_are_optimal(mesh));

    mesh.save_as(MDIR, "stator0");

    mesh.refine();
    mesh.save_as(MDIR, "stator0_refined");
}