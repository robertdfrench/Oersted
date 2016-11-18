#include "test_Sketch.hpp"

bool has_mirror_image(Sketch &s, std::shared_ptr<Vertex> v, LineSegment &l) {
    double x0 = l.start()->x();
    double y0 = l.start()->y();
    double x1 = l.end()->x();
    double y1 = l.end()->y();

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
        std::shared_ptr<Vertex> vj = s.vertex(j);
        mirror_found = mirror_found || (abs(vj->x() - vx) < TOL && abs(vj->y() - vy) < TOL);

        if (mirror_found) {
            break;
        }
    }

    return mirror_found;
}

void test_mirror_verticies(Sketch &s, std::vector<size_t> index, LineSegment &l) {
    for (size_t i : index) {
        EXPECT_TRUE(has_mirror_image(s, s.vertex(i), l));
    }
}

void test_mirror_curves(Sketch &s, std::vector<size_t> index, LineSegment &l) {
    for (size_t i : index) {
        const Curve *c = s.curve(i);
        std::shared_ptr<Vertex> v0 = c->start();
        std::shared_ptr<Vertex> v1 = c->end();

        bool mirror_found = has_mirror_image(s, v0, l) && has_mirror_image(s, v1, l);
        EXPECT_TRUE(mirror_found);
    }
}

TEST(MirrorCopy, nonoverlapping) {
    Sketch s;

    auto v0 = s.new_element_SHARED_PTR<Vertex>(2.0, 0.0);
    auto v1 = s.new_element_SHARED_PTR<Vertex>(3.0, 0.0);
    auto v2 = s.new_element_SHARED_PTR<Vertex>(3.0, 1.0);
    auto v3 = s.new_element_SHARED_PTR<Vertex>(2.0, 1.0);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
    LineSegment &l3 = s.new_element<LineSegment>(v3, v0);

    auto v4 = s.new_element_SHARED_PTR<Vertex>(-1.0, -1.0);
    auto v5 = s.new_element_SHARED_PTR<Vertex>(1.0, 1.0);
    LineSegment &l4 = s.new_element<LineSegment>(v4, v5);
    l4.ForConstruction = true;

    s.new_element<Fixation>(v4);
    s.new_element<Fixation>(v5);

    std::vector<const Curve *> vec{&l0, &l1, &l2, &l3};

    MirrorCopy &mc0 = s.new_element<MirrorCopy>(vec, &l4);

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Pattern__Mirror_nonoverlapping_square");

    s.solve();
    s.build();

    // Run Tests
    {
        test_sketch_size(s, 10, 9, 6, 2);
        test_mirror_verticies(s, {0, 1, 2, 3}, l4);
        test_mirror_curves(s, {0, 1, 2, 3}, l4);
    }

    // Change elements
    {
        s.new_element<Length>(l0, 0.5);
        s.solve();
        s.build();

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Pattern__Mirror_nonoverlapping_trapezoid");
    }

    // Run Tests
    {
        test_sketch_size(s, 10, 9, 7, 2);
        test_mirror_verticies(s, {0, 1, 2, 3}, l4);
        test_mirror_curves(s, {0, 1, 2, 3}, l4);
    }

    s.delete_me();
}

TEST(MirrorCopy, overlapping) {
    for (bool remove_internal : {true, false}) {
        Sketch s;

        auto v0 = s.new_element_SHARED_PTR<Vertex>(0.0, 0.0);
        auto v1 = s.new_element_SHARED_PTR<Vertex>(1.0, 0.0);
        auto v2 = s.new_element_SHARED_PTR<Vertex>(2.0, 1.0);
        auto v3 = s.new_element_SHARED_PTR<Vertex>(1.0, 1.0);

        LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
        LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
        LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
        LineSegment &l3 = s.new_element<LineSegment>(v3, v0);

        l3.ForConstruction = true;
        s.new_element<Fixation>(v0);
        s.new_element<Fixation>(v3);

        // std::vector<const Curve *> vec{&l0, &l1, &l2, &l3};
        auto mvec = s.curves();
        MirrorCopy &mc0 = s.new_element<MirrorCopy>(mvec, &l3, remove_internal);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Mirror__overlapping_parallelogram")+std::to_string(remove_internal));

        s.solve();
        s.build();

        // Run Tests
        {
            test_sketch_size(s, 6, 7, 4, 2 - remove_internal);
            test_mirror_verticies(s, {1, 2}, l3);
            test_mirror_curves(s, {0, 1, 2}, l3);
        }

        // Change elements
        {
            s.new_element<Length>(l1, 0.5);
            s.solve();
            s.build();

            s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Mirror__overlapping_trapezoid")+std::to_string(remove_internal));
        }

        // Run Tests
        {
            test_sketch_size(s, 6, 7, 5, 2 - remove_internal);
            test_mirror_verticies(s, {1, 2}, l3);
            test_mirror_curves(s, {0, 1, 2}, l3);
        }

        s.delete_me();
    }
}

TEST(MirrorCopy, multiple_overlapping) {
    for (bool remove_internal : {true, false}) {
        Sketch s;

        auto v0 = s.new_element_SHARED_PTR<Vertex>(0.0, 0.0);
        auto v1 = s.new_element_SHARED_PTR<Vertex>(2.0, 1.0);
        auto v2 = s.new_element_SHARED_PTR<Vertex>(2.0, 2.0);
        auto v3 = s.new_element_SHARED_PTR<Vertex>(1.0, 3.0);
        auto v4 = s.new_element_SHARED_PTR<Vertex>(3.0, 2.0);
        auto v5 = s.new_element_SHARED_PTR<Vertex>(2.0, 6.0);

        LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
        LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
        LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
        LineSegment &l3 = s.new_element<LineSegment>(v3, v0);
        LineSegment &l4 = s.new_element<LineSegment>(v2, v4);
        LineSegment &l5 = s.new_element<LineSegment>(v4, v5);
        LineSegment &l6 = s.new_element<LineSegment>(v5, v3);

        //l3.ForConstruction = true;
        s.new_element<Coincident<LineSegment>>(v5, l3);

        //std::vector<const Curve *> vec{&l0, &l1, &l2, &l3, &l4, &l5, &l6};
        auto mvec = s.curves();
        s.new_element<MirrorCopy>(mvec, &l3, remove_internal);

        s.solve();
        s.build();

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Mirror__multiple_overlapping_0_")+std::to_string(remove_internal));

        // Run Tests
        {
            test_sketch_size(s, 9, 12, 4, 4 - 2 * remove_internal);
            test_mirror_verticies(s, {1, 2, 4}, l3);
            test_mirror_curves(s, {0, 1, 2, 4, 5}, l3);
        }

        // Change elements
        {
            s.new_element<Length>(l6, 2.0);
            s.solve();
            s.build();

            s.save_as<SaveMethod::Rasterize>(SAVE_DIR, std::string("Mirror__multiple_overlapping_1_")+std::to_string(remove_internal));
        }

        // Run Tests
        {
            test_sketch_size(s, 9, 12, 5, 4 - 2 * remove_internal);
            test_mirror_verticies(s, {1, 2, 4}, l3);
            test_mirror_curves(s, {0, 1, 2, 4, 5}, l3);
        }

        s.delete_me();
    }
}