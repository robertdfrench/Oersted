#include "test_Sketch.hpp"

TEST(RotateCopy, nonoverlapping) {
    {
        Sketch s;

        size_t N = 4;
        double a_deg = 360.0 / N;
        double a_rad = M_PI * 2.0 / N;

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(1.0, 0.0);
        auto v2 = s.new_element<Vertex>(2.0, 0.0);
        auto v3 = s.new_element<Vertex>(2.0 * cos(a_rad), 2.0 * sin(a_rad));
        auto v4 = s.new_element<Vertex>(1.0 * cos(a_rad), 1.0 * sin(a_rad));

        auto l0 = s.new_element<LineSegment>(v1, v2);
        auto l1 = s.new_element<LineSegment>(v4, v3);

        auto c0 = s.new_element<CircularArc>(v2, v3, v0, 2.0);
        auto c1 = s.new_element<CircularArc>(v1, v4, v0, 1.0);

        auto rad0 = s.new_element<Radius>(c0, 2.0);
        auto rad1 = s.new_element<Radius>(c1, 1.0);

        auto f0 = s.new_element<Fixation>(v0);
        auto h0 = s.new_element<Horizontal>(l0);
        auto co0 = s.new_element<Coincident<LineSegment>>(v0, l1);
        auto a0 = s.new_element<Angle>(l0, l1, a_deg);

        std::vector<std::shared_ptr<Curve>> vec{l0, l1, c0, c1};
        auto r0 = s.new_element<RotateCopy>(vec, v0, 360.0 / (N - 1), N - 2);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Rotate__non_overlapping_0");

        s.solve();
        s.build();

        // Run Tests
        {
            test_sketch_size(s, 13, 12, 14, 3);
            test_rotated_verticies(s, {1, 2, 3, 4}, v0, 360.0 / (N - 1), N - 2);
            test_rotated_curves(s, {0, 1, 2, 3}, v0, 360.0 / (N - 1), N - 2);
        }

        // Change Sketch
        a0->Dim = a0->Dim / 2.0;
        rad0->Dim = 1.5;
        rad1->Dim = 0.5;

        s.solve();
        s.build();

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Rotate__non_overlapping_1");

        // Run Tests
        {
            test_sketch_size(s, 13, 12, 14, 3);
            test_rotated_verticies(s, {1, 2, 3, 4}, v0, 360.0 / (N - 1), N - 2);
            test_rotated_curves(s, {0, 1, 2, 3}, v0, 360.0 / (N - 1), N - 2);
        }
    }
}

TEST(RotateCopy, overlapping) {
    for (bool remove_internal : {false, true}) {
        Sketch s;

        size_t N = 4;
        double a_deg = 360.0 / N;
        double a_rad = M_PI * 2.0 / N;

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(1.0, 0.1);
        auto v2 = s.new_element<Vertex>(2.0, 0.2);
        auto v3 = s.new_element<Vertex>(2.0 * cos(a_rad), 2.0 * sin(a_rad));
        auto v4 = s.new_element<Vertex>(1.0 * cos(a_rad), 1.0 * sin(a_rad));
        auto v5 = s.new_element<Vertex>(3.0, 0.0);

        auto l0 = s.new_element<LineSegment>(v1, v2);
        auto l1 = s.new_element<LineSegment>(v4, v3);
        auto lh = s.new_element<LineSegment>(v0, v5);
        lh->ForConstruction = remove_internal;

        auto c0 = s.new_element<CircularArc>(v2, v3, v0, 2.0);
        auto c1 = s.new_element<CircularArc>(v1, v4, v0, 1.0);

        auto h = s.new_element<Horizontal>(lh);

        auto rad0 = s.new_element<Radius>(c0, 2.0);
        auto rad1 = s.new_element<Radius>(c1, 1.0);

        auto f0 = s.new_element<Fixation>(v0);

        auto co0 = s.new_element<Coincident<LineSegment>>(v0, l0);
        auto co1 = s.new_element<Coincident<LineSegment>>(v0, l1);

        auto a0 = s.new_element<Angle>(lh, l0, 22.5);
        //Horizontal &h = s.new_element<Horizontal>(l0);

        auto a1 = s.new_element<Angle>(l0, l1, a_deg);

        double res_norm = s.solve();
        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Rotate__overlapping_base_0_") +
                                                   std::to_string(remove_internal));

        bool result = s.build();
        ASSERT_TRUE(result);

        std::vector<std::shared_ptr<Curve>> rvec = {l0, l1, c0, c1};
        auto r0 = s.new_element<RotateCopy>(rvec, v0, 360.0 / N, N - 2, remove_internal);

        res_norm = s.solve();
        s.save_as<SaveMethod::Rasterize>(SAVE_DIR,
                                         std::string("Rotate__overlapping_0_") + std::to_string(remove_internal));

        result = s.build();
        ASSERT_TRUE(result);

        // Run Tests
        {
            test_sketch_size(s, 10, 11, 12, 3 - 2 * remove_internal);
            test_rotated_verticies(s, {1, 2, 3, 4}, v0, 360.0 / N, N - 2);
            test_rotated_curves(s, {0, 1, 3, 4}, v0, 360.0 / N, N - 2);
        }

        // Change Sketch
        rad0->Dim = 2.5;
        rad1->Dim = 0.5;
        a0->Dim = 45;

        res_norm = s.solve();
        s.save_as<SaveMethod::Rasterize>(SAVE_DIR,
                                         std::string("Rotate__overlapping_1_") + std::to_string(remove_internal));

        result = s.build();
        ASSERT_TRUE(result);

        // Run Tests
        {
            test_sketch_size(s, 10, 11, 12, 3 - 2 * remove_internal);
            test_rotated_verticies(s, {1, 2, 3, 4}, v0, 360.0 / N, N - 2);
            test_rotated_curves(s, {0, 1, 3, 4}, v0, 360.0 / N, N - 2);
        }
    }
}

TEST(RotateCopy, open_overlapping) {
    FAIL(); // TODO
}

TEST(RotateCopy, closed_overlapping) {
    FAIL(); // TODO: Need to implement detection of completely closed rotation copy
}

TEST(RotateCopy, noncoincident_overlapping) {
    FAIL(); // TODO: Current implementation will fail if rotated boundaries overlap but are not identical
}