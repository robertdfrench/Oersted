#include "test_Sketch.hpp"

bool has_mirror_image(Sketch &s, std::shared_ptr<Vertex const> v, std::shared_ptr<LineSegment> l) {
    double x0 = l->start()->x();
    double y0 = l->start()->y();
    double x1 = l->end()->x();
    double y1 = l->end()->y();

    double dx = x1 - x0;
    double dy = y1 - y0;
    double dl = sqrt(dx * dx + dy * dy);
    dx /= dl;
    dy /= dl;

    double dot = dx * (v->x() - x0) + dy * (v->y() - y0);
    double px = dot * dx;
    double py = dot * dy;
    double rx = (v->x() - x0) - px;
    double ry = (v->y() - y0) - py;

    double vx = px - rx + x0;
    double vy = py - ry + y0;

    bool mirror_found = false;
    for (size_t j = 0; j != s.size_verticies(); ++j) {
        std::shared_ptr<Vertex const> vj = s.vertex(j);
        mirror_found = mirror_found || (abs(vj->x() - vx) < TOL && abs(vj->y() - vy) < TOL);

        if (mirror_found) {
            break;
        }
    }

    return mirror_found;
}

void test_mirror_verticies(Sketch &s, std::vector<size_t> index, std::shared_ptr<LineSegment> l) {
    for (size_t i : index) {
        EXPECT_TRUE(has_mirror_image(s, s.vertex(i), l));
    }
}

void test_mirror_curves(Sketch &s, std::vector<size_t> index, std::shared_ptr<LineSegment> l) {
    for (size_t i : index) {
        auto c = s.curve(i);
        auto v0 = c->start();
        auto v1 = c->end();

        bool mirror_found = has_mirror_image(s, v0, l) && has_mirror_image(s, v1, l);
        EXPECT_TRUE(mirror_found);
    }
}

TEST(MirrorCopy, nonoverlapping) {
    Sketch s;

    auto v0 = s.new_element<Vertex>(2.0, 0.0);
    auto v1 = s.new_element<Vertex>(3.0, 0.0);
    auto v2 = s.new_element<Vertex>(3.0, 1.0);
    auto v3 = s.new_element<Vertex>(2.0, 1.0);

    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v1, v2);
    auto l2 = s.new_element<LineSegment>(v2, v3);
    auto l3 = s.new_element<LineSegment>(v3, v0);

    auto v4 = s.new_element<Vertex>(-1.0, -1.0);
    auto v5 = s.new_element<Vertex>(1.0, 1.0);
    auto l4 = s.new_element<LineSegment>(v4, v5);
    l4->for_construction(true);

    s.new_element<Fixation>(v4);
    s.new_element<Fixation>(v5);

    std::vector<std::shared_ptr<Curve const>> vec{l0, l1, l2, l3};

    auto mc0 = s.new_element<MirrorCopy>(vec, l4);

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Pattern__Mirror_nonoverlapping_square");

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    bool result = s.build();
    ASSERT_TRUE(result);

    // Run Tests
    test_sketch_size(s, 10, 9, 6, 2);
    test_mirror_verticies(s, {0, 1, 2, 3}, l4);
    test_mirror_curves(s, {0, 1, 2, 3}, l4);

    // Change elements
    s.new_element<Length>(l0, 0.5);

    res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    result = s.build();
    ASSERT_TRUE(result);

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Pattern__Mirror_nonoverlapping_trapezoid");

    // Run Tests
    test_sketch_size(s, 10, 9, 7, 2);
    test_mirror_verticies(s, {0, 1, 2, 3}, l4);
    test_mirror_curves(s, {0, 1, 2, 3}, l4);
}

TEST(MirrorCopy, overlapping) {
    for (bool remove_internal : {true, false}) {
        Sketch s;

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(1.0, 0.0);
        auto v2 = s.new_element<Vertex>(2.0, 1.0);
        auto v3 = s.new_element<Vertex>(1.0, 1.0);

        auto l0 = s.new_element<LineSegment>(v0, v1);
        auto l1 = s.new_element<LineSegment>(v1, v2);
        auto l2 = s.new_element<LineSegment>(v2, v3);
        auto l3 = s.new_element<LineSegment>(v3, v0);

        l3->for_construction(true);
        s.new_element<Fixation>(v0);
        s.new_element<Fixation>(v3);

        auto mvec = s.curves();
        auto mc0 = s.new_element<MirrorCopy>(mvec, l3, remove_internal);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Mirror__overlapping_parallelogram") + std::to_string(remove_internal));

        double res_norm = s.solve();
        EXPECT_LE(res_norm, FLT_EPSILON);

        bool result = s.build();
        ASSERT_TRUE(result);

        // Run Tests
        test_sketch_size(s, 6, 7, 4, 2 - remove_internal);
        test_mirror_verticies(s, {1, 2}, l3);
        test_mirror_curves(s, {0, 1, 2}, l3);

        // Change elements
        s.new_element<Length>(l1, 0.5);

        res_norm = s.solve();
        EXPECT_LE(res_norm, FLT_EPSILON);

        result = s.build();
        ASSERT_TRUE(result);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Mirror__overlapping_trapezoid") + std::to_string(remove_internal));

        // Run Tests
        test_sketch_size(s, 6, 7, 5, 2 - remove_internal);
        test_mirror_verticies(s, {1, 2}, l3);
        test_mirror_curves(s, {0, 1, 2}, l3);
    }
}

TEST(MirrorCopy, multiple_overlapping) {
    for (bool remove_internal : {true, false}) {
        Sketch s;

        auto v0 = s.new_element<Vertex>(0.0, 0.0);
        auto v1 = s.new_element<Vertex>(2.0, 1.0);
        auto v2 = s.new_element<Vertex>(2.0, 2.0);
        auto v3 = s.new_element<Vertex>(1.0, 3.0);
        auto v4 = s.new_element<Vertex>(3.0, 2.0);
        auto v5 = s.new_element<Vertex>(2.0, 6.0);

        auto l0 = s.new_element<LineSegment>(v0, v1);
        auto l1 = s.new_element<LineSegment>(v1, v2);
        auto l2 = s.new_element<LineSegment>(v2, v3);
        auto l3 = s.new_element<LineSegment>(v3, v0);
        auto l4 = s.new_element<LineSegment>(v2, v4);
        auto l5 = s.new_element<LineSegment>(v4, v5);
        auto l6 = s.new_element<LineSegment>(v5, v3);

        //l3.ForConstruction = true;
        s.new_element<Coincident<LineSegment>>(v5, l3);

        auto mvec = s.curves();
        s.new_element<MirrorCopy>(mvec, l3, remove_internal);

        double res_norm = s.solve();
        EXPECT_LE(res_norm, FLT_EPSILON);

        bool result = s.build();
        ASSERT_TRUE(result);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Mirror__multiple_overlapping_0_") + std::to_string(remove_internal));

        // Run Tests
        test_sketch_size(s, 9, 12, 4, 4 - 2 * remove_internal);
        test_mirror_verticies(s, {1, 2, 4}, l3);
        test_mirror_curves(s, {0, 1, 2, 4, 5}, l3);


        // Change elements
        s.new_element<Length>(l6, 2.0);

        res_norm = s.solve();
        EXPECT_LE(res_norm, FLT_EPSILON);

        result = s.build();
        ASSERT_TRUE(result);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Mirror__multiple_overlapping_1_") + std::to_string(remove_internal));

        // Run Tests
        test_sketch_size(s, 9, 12, 5, 4 - 2 * remove_internal);
        test_mirror_verticies(s, {1, 2, 4}, l3);
        test_mirror_curves(s, {0, 1, 2, 4, 5}, l3);
    }
}