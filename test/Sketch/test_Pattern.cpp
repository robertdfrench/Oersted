#include "test_Sketch.hpp"

void test_sketch_size(Sketch &s, size_t nverts, size_t ncurves, size_t nconstraints, size_t ncontours) {
    EXPECT_EQ(nverts, s.size_verticies());
    EXPECT_EQ(ncurves, s.size_curves());
    EXPECT_EQ(nconstraints, s.size_constraints());
    EXPECT_EQ(ncontours, s.size_contours());
}

bool has_mirror_image(Sketch &s, const Vertex *v, LineSegment &l) {
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
        const Vertex *vj = s.vertex(j);
        mirror_found = mirror_found || (abs(vj->x() - vx) < TOL && abs(vj->y() - vy) < TOL);

        if (mirror_found) {
            break;
        }
    }

    return mirror_found;
}

void test_mirror_verticies(Sketch &s, std::vector <size_t> index, LineSegment &l) {
    for (size_t i : index) {
        EXPECT_TRUE(has_mirror_image(s, s.vertex(i), l));
    }
}

void test_mirror_curves(Sketch &s, std::vector <size_t> index, LineSegment &l) {
    for (size_t i : index) {
        const Curve *c = s.curve(i);
        const Vertex *v0 = c->start();
        const Vertex *v1 = c->end();

        bool mirror_found = has_mirror_image(s, v0, l) && has_mirror_image(s, v1, l);
        EXPECT_TRUE(mirror_found);
    }
}

bool has_rotational_image(Sketch &s, const Vertex *v, const Vertex *center, double angle) {
    double x, y;
    std::tie(x, y) = v->rotate(center, angle);

    bool rotation_found = false;
    for (size_t j = 0; j != s.size_verticies(); ++j) {
        const Vertex *vj = s.vertex(j);
        rotation_found = (abs(vj->x() - x) < TOL && abs(vj->y() - y) < TOL);

        if (rotation_found) {
            break;
        }
    }

    return rotation_found;
}

bool has_rotational_image(Sketch &s, const Vertex *v0, const Vertex *v1, const Vertex *center, double angle) {
    double x0, y0;
    std::tie(x0, y0) = v0->rotate(center, angle);

    double x1, y1;
    std::tie(x1, y1) = v1->rotate(center, angle);

    bool rotation0_found = false;
    bool rotation1_found = false;

    for (size_t j = 0; j != s.size_curves(); ++j) {
        const Curve *c = s.curve(j);
        const Vertex *vs = c->start();
        const Vertex *ve = c->end();

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

void test_rotated_verticies(Sketch &s, std::vector <size_t> index, const Vertex *center, double angle, size_t copies) {
    for (size_t j = 1; j != (copies + 1); ++j) {
        for (size_t i : index) {
            EXPECT_TRUE(has_rotational_image(s, s.vertex(i), center, angle * j));
        }
    }
}

void test_rotated_curves(Sketch &s, std::vector <size_t> index, const Vertex *center, double angle, size_t copies) {
    for (size_t j = 1; j != (copies + 1); ++j) {
        for (size_t i : index) {
            EXPECT_TRUE(has_rotational_image(s, s.curve(i)->start(), s.curve(i)->end(), center, angle * j));
        }
    }
}

TEST(PATTERN, Mirror_nonoverlapping) {
        Sketch s;

        Vertex &v0 = s.new_element<Vertex>(2.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(3.0, 0.0);
        Vertex &v2 = s.new_element<Vertex>(3.0, 1.0);
        Vertex &v3 = s.new_element<Vertex>(2.0, 1.0);

        LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
        LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
        LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
        LineSegment &l3 = s.new_element<LineSegment>(v3, v0);

        Vertex &v4 = s.new_element<Vertex>(-1.0, -1.0);
        Vertex &v5 = s.new_element<Vertex>(1.0, 1.0);
        LineSegment &l4 = s.new_element<LineSegment>(v4, v5);
        l4.ForConstruction = true;

        s.new_element<Fixation>(v4);
        s.new_element<Fixation>(v5);

        std::vector<const Curve*> vec{ &l0, &l1, &l2, &l3 };

        MirrorCopy &mc0 = s.new_element<MirrorCopy>(vec, &l4);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "PATTERN__Mirror_nonoverlapping_square.csv");

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

            s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "PATTER__Mirror_nonoverlapping_trapezoid.csv");
        }

        // Run Tests
        {
            test_sketch_size(s, 10, 9, 7, 2);
            test_mirror_verticies(s, {0, 1, 2, 3}, l4);
            test_mirror_curves(s, {0, 1, 2, 3}, l4);
        }
}

TEST(PATTERN, MIRROR_overlapping) {
        Sketch s;

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(1.0, 0.0);
        Vertex &v2 = s.new_element<Vertex>(2.0, 1.0);
        Vertex &v3 = s.new_element<Vertex>(1.0, 1.0);

        LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
        LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
        LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
        LineSegment &l3 = s.new_element<LineSegment>(v3, v0);

        l3.ForConstruction = true;
        s.new_element<Fixation>(v0);
        s.new_element<Fixation>(v3);

        std::vector<const Curve*> vec{ &l0, &l1, &l2, &l3 };
        MirrorCopy &mc0 = s.new_element<MirrorCopy>(vec, &l3);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "PATTERN__Mirror_overlapping_parallelogram.csv");

        s.solve();
        s.build();

        // Run Tests
        {
            test_sketch_size(s, 6, 7, 4, 1);
            test_mirror_verticies(s, {1, 2}, l3);
            test_mirror_curves(s, {0, 1, 2}, l3);
        }

        // Change elements
        {
            s.new_element<Length>(l1, 0.5);
            s.solve();
            s.build();

            s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "PATTERN__Mirror_overlapping_trapezoid.csv");
        }

        // Run Tests
        {
            test_sketch_size(s, 6, 7, 5, 1);
            test_mirror_verticies(s, {1, 2}, l3);
            test_mirror_curves(s, {0, 1, 2}, l3);
        }
}

TEST(PATTERN, MIRROR_multiple_overlapping) {
        Sketch s;

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(2.0, 1.0);
        Vertex &v2 = s.new_element<Vertex>(2.0, 2.0);
        Vertex &v3 = s.new_element<Vertex>(1.0, 3.0);
        Vertex &v4 = s.new_element<Vertex>(3.0, 2.0);
        Vertex &v5 = s.new_element<Vertex>(2.0, 6.0);

        LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
        LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
        LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
        LineSegment &l3 = s.new_element<LineSegment>(v3, v0);
        LineSegment &l4 = s.new_element<LineSegment>(v2, v4);
        LineSegment &l5 = s.new_element<LineSegment>(v4, v5);
        LineSegment &l6 = s.new_element<LineSegment>(v5, v3);

        l3.ForConstruction = true;
        s.new_element<Coincident<LineSegment>>(v5, l3);

        std::vector<const Curve*> vec{ &l0, &l1, &l2, &l3, &l4, &l5, &l6 };
        s.new_element<MirrorCopy>(vec, &l3);

        s.solve();
        s.build();

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "PATTERN__Mirror_multiple_overlapping_0.csv");

        // Run Tests
        {
            test_sketch_size(s, 9, 12, 4, 3);
            test_mirror_verticies(s, {1, 2, 4}, l3);
            test_mirror_curves(s, {0, 1, 2, 4, 5}, l3);
        }

        // Change elements
        {
            s.new_element<Length>(l6, 2.0);
            s.solve();
            s.build();

            s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "PATTERN__Mirror_multiple_overlapping_1.csv");
        }

        // Run Tests
        {
            test_sketch_size(s, 9, 12, 5, 3);
            test_mirror_verticies(s, {1, 2, 4}, l3);
            test_mirror_curves(s, {0, 1, 2, 4, 5}, l3);
        }
}

TEST(PATTERN, Rotate_nonoverlapping) {
        Sketch s;

        size_t N = 4;
        double a_deg = 360.0 / N;
        double a_rad = M_PI * 2.0 / N;

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(1.0, 0.0);
        Vertex &v2 = s.new_element<Vertex>(2.0, 0.0);
        Vertex &v3 = s.new_element<Vertex>(2.0 * cos(a_rad), 2.0 * sin(a_rad));
        Vertex &v4 = s.new_element<Vertex>(1.0 * cos(a_rad), 1.0 * sin(a_rad));

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

        std::vector<const Curve*> vec{ &l0, &l1, &c0, &c1 };
        RotateCopy &r0 = s.new_element<RotateCopy>(vec, &v0, 360.0 / (N - 1), N - 2);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "PATTERN__Rotate_nonoverlapping_0.csv");

        s.solve();
        s.build();

        // Run Tests
        {
            test_sketch_size(s, 13, 12, 14, 3);
            test_rotated_verticies(s, {1, 2, 3, 4}, &v0, 360.0 / (N - 1), N - 2);
            test_rotated_curves(s, {0, 1, 2, 3}, &v0, 360.0 / (N - 1), N - 2);
        }

        // Change Sketch
        a0.Dim = a0.Dim / 2.0;
        rad0.Dim = 1.5;
        rad1.Dim = 0.5;

        s.solve();
        s.build();

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "PATTERN__Rotate_nonoverlapping_1.csv");

        // Run Tests
        {
            test_sketch_size(s, 13, 12, 14, 3);
            test_rotated_verticies(s, {1, 2, 3, 4}, &v0, 360.0 / (N - 1), N - 2);
            test_rotated_curves(s, {0, 1, 2, 3}, &v0, 360.0 / (N - 1), N - 2);
        }

        FAIL(); // #TODO: Rewrite rotate copy based on notes in current implementation
}

TEST(PATTERN, Rotate_overlapping) {
        Sketch s;

        size_t N = 4;
        double a_deg = 360.0 / N;
        double a_rad = M_PI * 2.0 / N;

        Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
        Vertex &v1 = s.new_element<Vertex>(1.0, 0.0);
        Vertex &v2 = s.new_element<Vertex>(2.0, 0.0);
        Vertex &v3 = s.new_element<Vertex>(2.0 * cos(a_rad), 2.0 * sin(a_rad));
        Vertex &v4 = s.new_element<Vertex>(1.0 * cos(a_rad), 1.0 * sin(a_rad));

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

        std::vector<const Curve*> vec({ &l0, &l1, &c0, &c1 });
        RotateCopy &r0 = s.new_element<RotateCopy>(vec, &v0, 360.0 / N, N - 2);

        s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "PATTERN__Rotate_overlapping_0.csv");

        s.solve();
        s.build();

        // Run Tests
        {
            test_sketch_size(s, 9, 8, 13, 1);
            test_rotated_verticies(s, {1, 2, 3, 4}, &v0, 360.0 / N, N - 2);
            test_rotated_curves(s, {0, 1, 2, 3}, &v0, 360.0 / N, N - 2);
        }

        // Change Sketch
        // Run Tests

        FAIL(); // #TODO: Change sketch, run tests
}

TEST(PATTERN, Rotate_overlapping_open) {
    FAIL(); // #TODO
}

TEST(PATTERN, Rotate_overlapping_closed) {
    FAIL(); // #TODO
}

TEST(PATTERN, Rotate_overlapping_noncoincident) {
    FAIL(); // #TODO: Current implementation will fail if rotated boundaries overlap but are not identical
}