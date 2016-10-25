#include "test_UseCases.hpp"

#define SAVE_DIR "./test/output/UseCases/Stator/"

TEST(STATOR, 0) {
    Sketch sketch;

    size_t Np{8};
    size_t Nt{3 * Np};

    double ri{80.0};
    double ro{130.0};
    double at{M_PI / Nt};
    double as{at / 2.0};
    double dt{ri * sin(as)};
    double rb{ro - 4.0 * dt};

    Vertex &origin = sketch.new_element<Vertex>(0.0, 0.0);
    Vertex &v0 = sketch.new_element<Vertex>(ri, 0.0);
    Vertex &v1 = sketch.new_element<Vertex>(ro, 0.0);
    Vertex &v2 = sketch.new_element<Vertex>(ro * cos(at), ro * sin(at));
    Vertex &v3 = sketch.new_element<Vertex>(rb * cos(at), rb * sin(at));
    Vertex &v4 = sketch.new_element<Vertex>(ri * cos(at), ri * sin(at));
    Vertex &v5 = sketch.new_element<Vertex>(ri * cos(as), ri * sin(as));
    Vertex &v6 = sketch.new_element<Vertex>(rb * cos(as), rb * sin(as));

    LineSegment &l0 = sketch.new_element<LineSegment>(v0, v1);
    CircularArc &c0 = sketch.new_element<CircularArc>(v1, v2, origin, ro);
    LineSegment &l1 = sketch.new_element<LineSegment>(v2, v3);
    LineSegment &l2 = sketch.new_element<LineSegment>(v3, v4);
    CircularArc &c1 = sketch.new_element<CircularArc>(v5, v4, origin, ri);
    CircularArc &c2 = sketch.new_element<CircularArc>(v0, v5, origin, ri);
    LineSegment &l3 = sketch.new_element<LineSegment>(v5, v6);
    CircularArc &c3 = sketch.new_element<CircularArc>(v6, v3, origin, rb);

    //std::vector<const Curve*> mv{ &l0,&l1,&l2,&l3,&c0,&c1,&c2,&c3 };
    std::vector<const Curve*> mv{&l0, &l3, &c0, &c1, &c2, &c3};
    // l1.ForConstruction = true; //TODO: Fixes segfault (#1) failure but should not be necessary, problem is related to MirrorCopy mirror line l1
    MirrorCopy &m0 = sketch.new_element<MirrorCopy>(mv, &l1);

    sketch.solve();

    sketch.save_as<SaveMethod::Rasterize>(SAVE_DIR, "stator0d0.csv");

    EXPECT_EQ(sketch.boundary()->size(), 8);

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

    sketch.save_as<SaveMethod::Rasterize>(SAVE_DIR, "stator0d1.csv");

    EXPECT_EQ(sketch.boundary()->size(), 8);

    sketch.build();

    Mesh mesh{sketch};

    mesh.MinimumElementQuality = M_SQRT1_2;
    mesh.MaximumElementSize = 2.5;
    mesh.MinimumElementSize = 0.25;

    /* TODO: Segfault (#1)
    mesh.create();

    mesh.save_as(SAVE_DIR, "stator_0.csv");

    mesh.refine();

    mesh.save_as(SAVE_DIR, "stator_0_refine.csv");
    */
}