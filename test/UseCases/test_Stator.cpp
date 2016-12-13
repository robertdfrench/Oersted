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

    auto origin = sketch.new_element<Vertex>(0.0, 0.0);
    auto v0 = sketch.new_element<Vertex>(ri, 0.0);
    auto v1 = sketch.new_element<Vertex>(ro, 0.0);
    auto v2 = sketch.new_element<Vertex>(ro * cos(at), ro * sin(at));
    auto v3 = sketch.new_element<Vertex>(rb * cos(at), rb * sin(at));
    auto v4 = sketch.new_element<Vertex>(ri * cos(at), ri * sin(at));
    auto v5 = sketch.new_element<Vertex>(ri * cos(as), ri * sin(as));
    auto v6 = sketch.new_element<Vertex>(rb * cos(as), rb * sin(as));

    auto l0 = sketch.new_element<LineSegment>(v0, v1);
    auto c0 = sketch.new_element<CircularArc>(v1, v2, origin, ro);
    auto l1 = sketch.new_element<LineSegment>(v2, v3);
    auto l2 = sketch.new_element<LineSegment>(v3, v4);
    auto c1 = sketch.new_element<CircularArc>(v5, v4, origin, ri);
    auto c2 = sketch.new_element<CircularArc>(v0, v5, origin, ri);
    auto l3 = sketch.new_element<LineSegment>(v5, v6);
    auto c3 = sketch.new_element<CircularArc>(v6, v3, origin, rb);

    auto mvec = sketch.curves();
    auto m0 = sketch.new_element<MirrorCopy>(mvec, l1, true);

    auto rvec = sketch.curves();
    auto rcopy = sketch.new_element<RotateCopy>(rvec, origin, 360.0 / Nt, 1, true);

    double res_norm = sketch.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    sketch.save_as<SaveMethod::Rasterize>(SDIR, "stator0_0");

    ASSERT_TRUE(sketch.build());
    EXPECT_EQ(sketch.boundary()->size(), 14);

    auto f = sketch.new_element<Fixation>(origin);

    auto coin0 = sketch.new_element<Coincident<LineSegment>>(origin, l0);
    auto coin1 = sketch.new_element<Coincident<LineSegment>>(origin, l1);
    auto coin2 = sketch.new_element<Coincident<LineSegment>>(origin, l2);

    auto horz0 = sketch.new_element<Horizontal>(l0);

    auto rad0 = sketch.new_element<Radius>(c0, ro);
    auto rad1 = sketch.new_element<Radius>(c1, ri);
    auto rad2 = sketch.new_element<Radius>(c2, ri);
    auto rad3 = sketch.new_element<Radius>(c3, rb);

    auto ang0 = sketch.new_element<Angle>(l1, l0, (at * 180.0 / M_PI));
    auto ang1 = sketch.new_element<Angle>(l2, l0, (at * 180.0 / M_PI));

    auto dist0 = sketch.new_element<Distance<LineSegment>>(l0, l3, dt);

    res_norm = sketch.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    sketch.save_as<SaveMethod::Rasterize>(SDIR, "stator0_1");

    ASSERT_TRUE(sketch.build());
    EXPECT_EQ(sketch.boundary()->size(), 14);

    Mesh mesh{sketch};

    mesh.MinimumElementQuality = 0.1;
    mesh.MaximumElementSize = 4.0;
    mesh.MinimumElementSize = 0.40;

    mesh.create();

    EXPECT_TRUE(edges_are_valid(mesh));
    EXPECT_TRUE(edges_are_optimal(mesh));

    mesh.save_as(MDIR, "stator0");

    mesh.refine();
    mesh.save_as(MDIR, "stator0_refined");
}

TEST(Stator, 1) {
    Sketch sketch;

    // Parameters
    size_t Np{8}; // number of poles
    size_t Nt{6 * Np}; // number of teeth

    double rsi{80e-3}; // stator inner radius
    double rso{130e-3}; // stator outer radius
    double at2_deg{180.0 / Nt}; // half tooth pitch angle in degrees
    double at2_rad{M_PI / Nt}; // half tooth pitch angle in radians

    double sow{1.0e-3}; // slot opening width
    double sod{1.5e-3}; // slot opening depth

    double ptw{0.66}; // tooth width percentage of rsi * cos(at2_rad)
    double dtw{ptw * rsi * sin(at2_rad)};

    double pbi{0.3}; // backiron width percentage of (rso - rsi)
    double dbi{pbi * (rso - rsi)};

    // Verticies
    auto v0 = sketch.new_element<Vertex>(0.0, 0.0);
    auto f0 = sketch.new_element<Fixation>(v0);

    auto v1 = sketch.new_element<Vertex>(rsi, 0.0);
    auto v2 = sketch.new_element<Vertex>(rso, 0.0);
    auto v3 = sketch.new_element<Vertex>(rso * cos(at2_rad), rso * sin(at2_rad));
    auto v4 = sketch.new_element<Vertex>(rsi * cos(at2_rad), rsi * sin(at2_rad));

    auto v5 = sketch.new_element<Vertex>(rsi * cos(at2_rad), rsi * sin(at2_rad) - sod); // approximate dimensions
    auto v6 = sketch.new_element<Vertex>(rsi * cos(at2_rad) + sod, rsi * sin(at2_rad) - sod); // approximate dimensions
    auto v7 = sketch.new_element<Vertex>((rsi + sod) * cos(at2_rad), (rsi + sod) * sin(at2_rad));

    auto v8 = sketch.new_element<Vertex>(rsi + 2 * sod, dtw);
    auto v9 = sketch.new_element<Vertex>(rso - dbi, dtw);
    auto v10 = sketch.new_element<Vertex>((rso - dbi) * cos(at2_rad), (rso - dbi) * sin(at2_rad));

    auto v11 = sketch.new_element<Vertex>(rsi + 2 * sod, rsi * sin(at2_rad) - sod);
    auto v12 = sketch.new_element<Vertex>((rso - 2 * dbi) * cos(at2_rad), (rso - 2 * dbi) * sin(at2_rad));

    // LineSegment Curves
    auto l_0_4 = sketch.new_element<LineSegment>(v0, v4, true);

    auto l_1_2 = sketch.new_element<LineSegment>(v1, v2);
    auto l_4_7 = sketch.new_element<LineSegment>(v4, v7);
    auto l_7_10 = sketch.new_element<LineSegment>(v7, v10);
    auto l_10_3 = sketch.new_element<LineSegment>(v10, v3);

    auto l_5_6 = sketch.new_element<LineSegment>(v5, v6);
    auto l_6_7 = sketch.new_element<LineSegment>(v6, v7);

    auto l_8_9 = sketch.new_element<LineSegment>(v8, v9);

    // Circle Curves
    auto c_1_5_0 = sketch.new_element<CircularArc>(v1, v5, v0, rsi);
    auto c_5_4_0 = sketch.new_element<CircularArc>(v5, v4, v0, rsi);
    auto c_2_3_0 = sketch.new_element<CircularArc>(v2, v3, v0, rso);

    auto c_6_8_11 = sketch.new_element<CircularArc>(v6, v8, v11);

    auto c_9_10_12 = sketch.new_element<CircularArc>(v9, v10, v12);

    // Angle Constraints
    auto angle_1_2_0_4 = sketch.new_element<Angle>(l_1_2, l_0_4, at2_deg);
    auto angle_1_2_4_7 = sketch.new_element<Angle>(l_1_2, l_4_7, at2_deg);
    auto angle_1_2_7_10 = sketch.new_element<Angle>(l_1_2, l_7_10, at2_deg);
    auto angle_1_2_10_3 = sketch.new_element<Angle>(l_1_2, l_10_3, at2_deg);

    // Distance Constraint
    auto distance_5_6_4_7 = sketch.new_element<Distance<LineSegment>>(l_5_6, l_4_7, sow);
    auto distance_1_2_8_9 = sketch.new_element<Distance<LineSegment>>(l_1_2, l_8_9, dtw);

    // Coincident Constraints
    auto coincident_0_1_2 = sketch.new_element<Coincident<LineSegment>>(v0, l_1_2);

    auto coincident_12_7_10 = sketch.new_element<Coincident<LineSegment>>(v12, l_7_10);

    // Horizontal Constraints
    auto horziontal_1_2 = sketch.new_element<Horizontal>(l_1_2);

    // Length Constraints
    auto length_5_6 = sketch.new_element<Length>(l_5_6, sod);
    auto length_4_7 = sketch.new_element<Length>(l_4_7, sod);
    auto length_10_3 = sketch.new_element<Length>(l_10_3, dbi);

    // Radius Constraints
    auto radius_1_5_0 = sketch.new_element<Radius>(c_1_5_0, rsi);
    auto radius_2_3_0 = sketch.new_element<Radius>(c_2_3_0, rso);

    // Tangency Constraints
    auto tangent_9_10_12_8_9 = sketch.new_element<Tangency>(c_9_10_12, l_8_9);

    auto tangent_6_8_11_8_9 = sketch.new_element<Tangency>(c_6_8_11, l_8_9);

    auto tangent_6_8_11_6_7 = sketch.new_element<Tangency>(c_6_8_11, l_6_7);

    // Solve
    double res_norm = sketch.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    sketch.save_as<SaveMethod::Rasterize>(SDIR, "stator1_0_half");

    // MirrorCopy
    auto mccs = sketch.curves();

    auto mirror_copy = sketch.new_element<MirrorCopy>(mccs, l_0_4, true);

    // Solve Again
    res_norm = sketch.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    sketch.save_as<SaveMethod::Rasterize>(SDIR, "stator1_1_mirror");

    // RotateCopy
    auto rccs = sketch.curves();

    auto rotate_copy = sketch.new_element<RotateCopy>(rccs, v0, at2_deg * 2.0, 1, true); // Can do 5 copies for full pole, but is slow at -O0

    // Solve Again
    res_norm = sketch.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    sketch.save_as<SaveMethod::Rasterize>(SDIR, "stator1_2_rotate");

    // Solve Again
    radius_1_5_0->dim(radius_1_5_0->dim() + 10e-3);

    res_norm = sketch.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    sketch.save_as<SaveMethod::Rasterize>(SDIR, "stator1_2_rotate_radius");

    // Build
    bool result = sketch.build();
    ASSERT_TRUE(result);

    // Create Mesh
    Mesh mesh{sketch};

    mesh.MinimumElementQuality = 0.1;
    mesh.MaximumElementSize = M_PI * radius_1_5_0->dim() / Nt / 3.0;
    mesh.MinimumElementSize = mesh.MaximumElementSize / 10.0;

    mesh.create();

    EXPECT_TRUE(edges_are_valid(mesh));
    EXPECT_TRUE(edges_are_optimal(mesh));

    mesh.save_as(MDIR, "stator1");

    mesh.refine();
    mesh.save_as(MDIR, "stator1_refined");
}