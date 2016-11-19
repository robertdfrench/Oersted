#include "test_Sketch.hpp"

bool has_rotational_image(Sketch &s, std::shared_ptr<Vertex> v, std::shared_ptr<Vertex> center, double angle) {
    double x, y;
    std::tie(x, y) = v->rotate(center, angle);

    bool rotation_found = false;
    for (size_t j = 0; j != s.size_verticies(); ++j) {
        std::shared_ptr<Vertex> vj = s.vertex(j);
        rotation_found = (abs(vj->x() - x) < TOL && abs(vj->y() - y) < TOL);

        if (rotation_found) {
            break;
        }
    }

    return rotation_found;
}

bool has_rotational_image(Sketch &s, std::shared_ptr<Vertex> v0, std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> center, double angle) {
    double x0, y0;
    std::tie(x0, y0) = v0->rotate(center, angle);

    double x1, y1;
    std::tie(x1, y1) = v1->rotate(center, angle);

    bool rotation0_found = false;
    bool rotation1_found = false;

    for (size_t j = 0; j != s.size_curves(); ++j) {
        const Curve *c = s.curve(j);
        std::shared_ptr<Vertex> vs = c->start();
        std::shared_ptr<Vertex> ve = c->end();

        rotation0_found = (abs(vs->x() - x0) < TOL && abs(vs->y() - y0) < TOL);

        if (rotation0_found) {
            rotation1_found = (abs(ve->x() - x1) < TOL && abs(ve->y() - y1) < TOL);
        } else {
            rotation0_found = (abs(ve->x() - x0) < TOL && abs(ve->y() - y0) < TOL);
            if (rotation0_found) {
                rotation1_found = (abs(vs->x() - x1) < TOL && abs(vs->y() - y1) < TOL);
            }
        }

        if (rotation0_found && rotation1_found) {
            break;
        }
    }

    return rotation0_found && rotation1_found;
}

void test_rotated_verticies(Sketch &s, std::vector<size_t> index, std::shared_ptr<Vertex> center, double angle, size_t copies) {
    for (size_t j = 1; j != (copies + 1); ++j) {
        for (size_t i : index) {
            EXPECT_TRUE(has_rotational_image(s, s.vertex(i), center, angle * j));
        }
    }
}

void test_rotated_curves(Sketch &s, std::vector<size_t> index, std::shared_ptr<Vertex> center, double angle, size_t copies) {
    for (size_t j = 1; j != (copies + 1); ++j) {
        for (size_t i : index) {
            EXPECT_TRUE(has_rotational_image(s, s.curve(i)->start(), s.curve(i)->end(), center, angle * j));
        }
    }
}

TEST(RotateCopy, nonoverlapping) {
    Sketch s;

    size_t N = 4;
    double a_deg = 360.0 / N;
    double a_rad = M_PI * 2.0 / N;

    auto v0 = s.new_element_SHARED_PTR<Vertex>(0.0, 0.0);
    auto v1 = s.new_element_SHARED_PTR<Vertex>(1.0, 0.0);
    auto v2 = s.new_element_SHARED_PTR<Vertex>(2.0, 0.0);
    auto v3 = s.new_element_SHARED_PTR<Vertex>(2.0 * cos(a_rad), 2.0 * sin(a_rad));
    auto v4 = s.new_element_SHARED_PTR<Vertex>(1.0 * cos(a_rad), 1.0 * sin(a_rad));

    LineSegment &l0 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l1 = s.new_element<LineSegment>(v4, v3);

    CircularArc &c0 = s.new_element<CircularArc>(v2, v3, v0, 2.0);
    CircularArc &c1 = s.new_element<CircularArc>(v1, v4, v0, 1.0);

    Radius &rad0 = s.new_element<Radius>(c0, 2.0);
    Radius &rad1 = s.new_element<Radius>(c1, 1.0);

    Fixation &f0 = s.new_element<Fixation>(v0);
    Horizontal &h0 = s.new_element<Horizontal>(l0);
    Coincident<LineSegment> &co0 = s.new_element<Coincident<LineSegment>>(v0, l1);
    Angle &a0 = s.new_element<Angle>(l0, l1, a_deg);

    std::vector<const Curve *> vec{&l0, &l1, &c0, &c1};
    auto r0 = s.new_element_SHARED_PTR<RotateCopy>(vec, v0, 360.0 / (N - 1), N - 2);

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
    a0.Dim = a0.Dim / 2.0;
    rad0.Dim = 1.5;
    rad1.Dim = 0.5;

    s.solve();
    s.build();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Rotate__non_overlapping_1");

    // Run Tests
    {
        test_sketch_size(s, 13, 12, 14, 3);
        test_rotated_verticies(s, {1, 2, 3, 4}, v0, 360.0 / (N - 1), N - 2);
        test_rotated_curves(s, {0, 1, 2, 3}, v0, 360.0 / (N - 1), N - 2);
    }

    s.delete_me();
}

TEST(RotateCopy, overlapping) {
    for (bool remove_internal : {true, false}) {
        Sketch s = Sketch();

        size_t N = 4;
        double a_deg = 360.0 / N;
        double a_rad = M_PI * 2.0 / N;

        auto v0 = s.new_element_SHARED_PTR<Vertex>(0.0, 0.0);
        auto v1 = s.new_element_SHARED_PTR<Vertex>(1.0, 0.1);
        auto v2 = s.new_element_SHARED_PTR<Vertex>(2.0, 0.2);
        auto v3 = s.new_element_SHARED_PTR<Vertex>(2.0 * cos(a_rad), 2.0 * sin(a_rad));
        auto v4 = s.new_element_SHARED_PTR<Vertex>(1.0 * cos(a_rad), 1.0 * sin(a_rad));
        auto v5 = s.new_element_SHARED_PTR<Vertex>(3.0, 0.0);

        LineSegment &l0 = s.new_element<LineSegment>(v1, v2);
        LineSegment &l1 = s.new_element<LineSegment>(v4, v3);
        LineSegment &lh = s.new_element<LineSegment>(v0, v5);
        lh.ForConstruction = remove_internal;

        CircularArc &c0 = s.new_element<CircularArc>(v2, v3, v0, 2.0);
        CircularArc &c1 = s.new_element<CircularArc>(v1, v4, v0, 1.0);

        Horizontal &h = s.new_element<Horizontal>(lh);

        Radius &rad0 = s.new_element<Radius>(c0, 2.0);
        Radius &rad1 = s.new_element<Radius>(c1, 1.0);

        Fixation &f0 = s.new_element<Fixation>(v0);

        Coincident<LineSegment> &co0 = s.new_element<Coincident<LineSegment>>(v0, l0);
        Coincident<LineSegment> &co1 = s.new_element<Coincident<LineSegment>>(v0, l1);

        Angle &a0 = s.new_element<Angle>(lh, l0, 22.5);
        //Horizontal &h = s.new_element<Horizontal>(l0);

        Angle &a1 = s.new_element<Angle>(l0, l1, a_deg);

        s.solve();

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Rotate__overlapping_base_0_") + std::to_string(remove_internal));

        std::vector<const Curve *> rvec = {&l0, &l1, &c0, &c1};
        auto r0 = s.new_element_SHARED_PTR<RotateCopy>(rvec, v0, 360.0 / N, N - 2, remove_internal);

        s.solve();
        s.build();

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Rotate__overlapping_0_") + std::to_string(remove_internal));

        // Run Tests
        {
            test_sketch_size(s, 10, 11, 12, 3 - 2 * remove_internal);
            test_rotated_verticies(s, {1, 2, 3, 4}, v0, 360.0 / N, N - 2);
            test_rotated_curves(s, {0, 1, 3, 4}, v0, 360.0 / N, N - 2);
        }

        // Change Sketch
        rad0.Dim = 2.5;
        rad1.Dim = 0.5;
        a0.Dim = 45;

        s.solve();
        s.build();

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Rotate__overlapping_1_") + std::to_string(remove_internal));

        // Run Tests
        {
            test_sketch_size(s, 10, 11, 12, 3 - 2 * remove_internal);
            test_rotated_verticies(s, {1, 2, 3, 4}, v0, 360.0 / N, N - 2);
            test_rotated_curves(s, {0, 1, 3, 4}, v0, 360.0 / N, N - 2);
        }

        s.delete_me();
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