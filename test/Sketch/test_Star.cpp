#include "test_Sketch.hpp"

bool star_angle_sum_equals_2pi(vector<Star> &stars) {
    for (auto s = stars.begin(); s != stars.end(); ++s) {
        double angle = 0.0;
        for (auto b = s->begin(); b != s->end(); ++b) {
            angle += b->Angle;
        }
        EXPECT_NEAR(2.0 * M_PI, angle, TOL);
    }

    return true;
}

TEST(Star, Suite_0) {
    Sketch sketch;
    Vertex &v0 = sketch.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = sketch.new_element<Vertex>(1.0, 1.0);
    Vertex &v2 = sketch.new_element<Vertex>(1.0, -1.0);
    Vertex &v3 = sketch.new_element<Vertex>(-1.0, 1.0);
    Vertex &v4 = sketch.new_element<Vertex>(-1.0, -1.0);

    LineSegment &l0 = sketch.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = sketch.new_element<LineSegment>(v2, v0);
    LineSegment &l2 = sketch.new_element<LineSegment>(v3, v0);
    LineSegment &l3 = sketch.new_element<LineSegment>(v0, v4);
    LineSegment &l4 = sketch.new_element<LineSegment>(v1, v2);
    LineSegment &l5 = sketch.new_element<LineSegment>(v2, v3);
    LineSegment &l6 = sketch.new_element<LineSegment>(v3, v4);
    LineSegment &l7 = sketch.new_element<LineSegment>(v4, v1);

    Star star{&v0, &sketch};

    EXPECT_TRUE(star.vertex() == &v0);
    EXPECT_TRUE(star.size() == 4);

    auto b = star.begin();

    EXPECT_TRUE((b->Path == &l2));
    EXPECT_TRUE(b->Orientation == false);

    EXPECT_TRUE((++b)->Path == &l0);
    EXPECT_TRUE(b->Orientation == true);

    EXPECT_TRUE((++b)->Path == &l1);
    EXPECT_TRUE(b->Orientation == false);

    EXPECT_TRUE((++b)->Path == &l3);
    EXPECT_TRUE(b->Orientation == true);

    for (auto b = star.begin(); b != star.end(); ++b) {
        //size_t j = (i + 1) % star.size();
        auto c = star.next(b);
        EXPECT_TRUE(star.next(b->Path) == c->Path);
    }

    EXPECT_TRUE(star.next(&l4) == nullptr);
    EXPECT_TRUE(star.next(&l5) == nullptr);
    EXPECT_TRUE(star.next(&l6) == nullptr);
    EXPECT_TRUE(star.next(&l7) == nullptr);

    double angle = 0.0;
    for (auto b = star.begin(); b != star.end(); ++b) {
        angle += b->Angle;
    }
    EXPECT_NEAR(2.0 * M_PI, angle, TOL);
}

TEST(Star, Suite_1) {
    Sketch sketch;

    Vertex &vs = sketch.new_element<Vertex>(0.0, 1.0);
    Vertex &vc = sketch.new_element<Vertex>(0.0, 0.0);
    Vertex &v0 = sketch.new_element<Vertex>(1.0, 0.0);
    Vertex &v1 = sketch.new_element<Vertex>(0.0, 2.0);
    Vertex &v2 = sketch.new_element<Vertex>(-1.0, 0.0);

    CircularArc &c0 = sketch.new_element<CircularArc>(v0, vs, vc, 1.0);
    CircularArc &c1 = sketch.new_element<CircularArc>(vs, v2, vc, 1.0);
    LineSegment &l0 = sketch.new_element<LineSegment>(vs, v1);
    LineSegment &l1 = sketch.new_element<LineSegment>(vc, vs);

    Star star{&vs, &sketch};

    EXPECT_TRUE(star.vertex() == &vs);
    EXPECT_TRUE(star.size() == 4);

    auto b = star.begin();

    EXPECT_TRUE(b->Path == &l0);
    EXPECT_TRUE(b->Orientation == true);

    EXPECT_TRUE((++b)->Path == &c0);
    EXPECT_TRUE(b->Orientation == false);

    EXPECT_TRUE((++b)->Path == &l1);
    EXPECT_TRUE(b->Orientation == false);

    EXPECT_TRUE((++b)->Path == &c1);
    EXPECT_TRUE(b->Orientation == true);

    for (auto b = star.begin(); b != star.end(); ++b) {
        //size_t j = (i + 1) % star.size();
        auto c = star.next(b);
        EXPECT_TRUE(star.next(b->Path) == c->Path);
    }

    EXPECT_TRUE(star.next(nullptr) == nullptr);

    double angle = 0.0;
    for (auto b = star.begin(); b != star.end(); ++b) {
        angle += b->Angle;
    }
    EXPECT_NEAR(2.0 * M_PI, angle, TOL);
}

TEST(Star, Suite_2) {
    Sketch sketch;

    Vertex &vs = sketch.new_element<Vertex>(M_SQRT1_2, M_SQRT1_2);
    Vertex &vc = sketch.new_element<Vertex>(0.0, 0.0);
    Vertex &v0 = sketch.new_element<Vertex>(1.0, 0.0);
    Vertex &v1 = sketch.new_element<Vertex>(0.0, 1.0);
    Vertex &v2 = sketch.new_element<Vertex>(M_SQRT2, 0.0);
    Vertex &v3 = sketch.new_element<Vertex>(0.0, M_SQRT2);
    Vertex &v4 = sketch.new_element<Vertex>(M_SQRT2, M_SQRT2);

    CircularArc &c0 = sketch.new_element<CircularArc>(v0, vs, vc, 1.0);
    CircularArc &c1 = sketch.new_element<CircularArc>(vs, v1, vc, 1.0);
    LineSegment &l0 = sketch.new_element<LineSegment>(vs, v2);
    LineSegment &l1 = sketch.new_element<LineSegment>(v3, vs);
    LineSegment &l2 = sketch.new_element<LineSegment>(vc, vs);
    LineSegment &l3 = sketch.new_element<LineSegment>(vs, v4);

    Star star{&vs, &sketch};

    EXPECT_TRUE(star.vertex() == &vs);
    EXPECT_TRUE(star.size() == 6);

    auto b = star.begin();

    EXPECT_TRUE(b->Path == &c1);
    EXPECT_TRUE(b->Orientation == true);

    EXPECT_TRUE((++b)->Path == &l1);
    EXPECT_TRUE(b->Orientation == false);

    EXPECT_TRUE((++b)->Path == &l3);
    EXPECT_TRUE(b->Orientation == true);

    EXPECT_TRUE((++b)->Path == &l0);
    EXPECT_TRUE(b->Orientation == true);

    EXPECT_TRUE((++b)->Path == &c0);
    EXPECT_TRUE(b->Orientation == false);

    EXPECT_TRUE((++b)->Path == &l2);
    EXPECT_TRUE(b->Orientation == false);

    for (auto b = star.begin(); b != star.end(); ++b) {
        //size_t j = (i + 1) % star.size();
        auto c = star.next(b);
        EXPECT_TRUE(star.next(b->Path) == c->Path);
    }

    EXPECT_TRUE(star.next(nullptr) == nullptr);

    double angle = 0.0;
    for (auto b = star.begin(); b != star.end(); ++b) {
        angle += b->Angle;
    }
    EXPECT_NEAR(2.0 * M_PI, angle, TOL);
}

TEST(Star, find_closed_contour_0) {
    Sketch sketch;
    Vertex &v0 = sketch.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = sketch.new_element<Vertex>(1.0, 0.0);
    Vertex &v2 = sketch.new_element<Vertex>(2.0, 0.0);
    Vertex &v3 = sketch.new_element<Vertex>(1.0, -1.0);
    Vertex &v4 = sketch.new_element<Vertex>(1.0, 1.0);
    Vertex &v5 = sketch.new_element<Vertex>(2.0, 1.0);
    Vertex &v6 = sketch.new_element<Vertex>(2.0, 2.0);
    Vertex &v7 = sketch.new_element<Vertex>(0.0, -1.0);
    Vertex &v8 = sketch.new_element<Vertex>(-1.0, 0.0);

    LineSegment &l0 = sketch.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = sketch.new_element<LineSegment>(v0, v7);
    LineSegment &l2 = sketch.new_element<LineSegment>(v0, v8);

    LineSegment &l3 = sketch.new_element<LineSegment>(v1, v2);
    LineSegment &l4 = sketch.new_element<LineSegment>(v1, v3);
    LineSegment &l5 = sketch.new_element<LineSegment>(v1, v4);

    LineSegment &l6 = sketch.new_element<LineSegment>(v4, v5);
    LineSegment &l7 = sketch.new_element<LineSegment>(v4, v6);
    LineSegment &l8 = sketch.new_element<LineSegment>(v4, v0);

    // Manual contour creation
    {
        vector<Star> stars;

        stars.push_back(Star{&v0, &sketch});
        stars.push_back(Star{&v1, &sketch});
        stars.push_back(Star{&v4, &sketch});

        star_angle_sum_equals_2pi(stars);
    }

    // Sketch internal contour parsing
    {
        sketch.solve();
        EXPECT_TRUE(sketch.build());
        EXPECT_EQ(sketch.size_contours(), 1);
        EXPECT_TRUE(*sketch.contour(0) == *sketch.boundary());

        const Contour *contour = sketch.contour(0);
        EXPECT_TRUE(contour->size() == 3);
        EXPECT_TRUE(&l0 == contour->curve(0) || &l0 == contour->curve(1) || &l0 == contour->curve(2));
        EXPECT_TRUE(&l5 == contour->curve(0) || &l5 == contour->curve(1) || &l5 == contour->curve(2));
        EXPECT_TRUE(&l8 == contour->curve(0) || &l8 == contour->curve(1) || &l8 == contour->curve(2));
    }
}

TEST(Star, find_closed_contour_1) {
    Sketch sketch;

    Vertex &v0 = sketch.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = sketch.new_element<Vertex>(1.0, 0.0);
    Vertex &v2 = sketch.new_element<Vertex>(0.0, 1.0);
    Vertex &v3 = sketch.new_element<Vertex>(0.0, -1.0);
    Vertex &v4 = sketch.new_element<Vertex>(-1.0, 0.0);

    LineSegment &l0 = sketch.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = sketch.new_element<LineSegment>(v0, v2);
    LineSegment &l2 = sketch.new_element<LineSegment>(v4, v0);

    CircularArc &c0 = sketch.new_element<CircularArc>(v1, v2, v0, 1.0);
    CircularArc &c1 = sketch.new_element<CircularArc>(v3, v1, v0, 1.0);

    // Manual contour construction
    {
        vector<Star> stars;
        stars.push_back(Star{&v0, &sketch});
        stars.push_back(Star{&v1, &sketch});
        stars.push_back(Star{&v2, &sketch});

        star_angle_sum_equals_2pi(stars);
    }

    // Sketch internal contour parsing
    {
        sketch.solve();
        sketch.build();
        EXPECT_TRUE(sketch.size_contours() == 1);
        EXPECT_TRUE(*sketch.contour(0) == *sketch.boundary());

        const Contour *contour = sketch.contour(0);
        EXPECT_TRUE(contour->size() == 3);
        EXPECT_TRUE(&l0 == contour->curve(0) || &l0 == contour->curve(1) || &l0 == contour->curve(2));
        EXPECT_TRUE(&l1 == contour->curve(0) || &l1 == contour->curve(1) || &l1 == contour->curve(2));
        EXPECT_TRUE(&c0 == contour->curve(0) || &c0 == contour->curve(1) || &c0 == contour->curve(2));
    }
}

TEST(Star, find_closed_contour_2) {
    Sketch sketch;

    Vertex &v0 = sketch.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = sketch.new_element<Vertex>(1.0, 0.0);
    Vertex &v2 = sketch.new_element<Vertex>(0.0, 1.0);
    Vertex &v3 = sketch.new_element<Vertex>(0.0, -1.0);
    Vertex &v4 = sketch.new_element<Vertex>(-1.0, 0.0);
    Vertex &v5 = sketch.new_element<Vertex>(1.0, 1.0);

    LineSegment &l0 = sketch.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = sketch.new_element<LineSegment>(v0, v2);
    LineSegment &l2 = sketch.new_element<LineSegment>(v4, v0);
    LineSegment &l3 = sketch.new_element<LineSegment>(v1, v5);
    LineSegment &l4 = sketch.new_element<LineSegment>(v5, v2);

    CircularArc &c0 = sketch.new_element<CircularArc>(v1, v2, v0, 1.0);
    CircularArc &c1 = sketch.new_element<CircularArc>(v3, v1, v0, 1.0);

    // Manual contour construction
    {
        vector<Star> stars;
        stars.push_back(Star{&v0, &sketch});
        stars.push_back(Star{&v1, &sketch});
        stars.push_back(Star{&v2, &sketch});
        stars.push_back(Star{&v5, &sketch});

        star_angle_sum_equals_2pi(stars);
    }

    // Sketch internal contour parsing
    {
        sketch.solve();
        EXPECT_TRUE(sketch.build());

        EXPECT_EQ(sketch.size_contours(), 2);
        EXPECT_FALSE(sketch.contour(0) == sketch.boundary());
        EXPECT_FALSE(sketch.contour(1) == sketch.boundary());

        const Contour *contour = sketch.contour(1);

        EXPECT_EQ(contour->size(), 3);
        EXPECT_TRUE(&l0 == contour->curve(0) || &l0 == contour->curve(1) || &l0 == contour->curve(2));
        EXPECT_TRUE(&l1 == contour->curve(0) || &l1 == contour->curve(1) || &l1 == contour->curve(2));
        EXPECT_TRUE(&c0 == contour->curve(0) || &c0 == contour->curve(1) || &c0 == contour->curve(2));

        contour = sketch.contour(0);

        EXPECT_EQ(contour->size(), 3);
        EXPECT_TRUE(&l3 == contour->curve(0) || &l3 == contour->curve(1) || &l3 == contour->curve(2));
        EXPECT_TRUE(&l4 == contour->curve(0) || &l4 == contour->curve(1) || &l4 == contour->curve(2));
        EXPECT_TRUE(&c0 == contour->curve(0) || &c0 == contour->curve(1) || &c0 == contour->curve(2));

        const Contour *c0 = sketch.contour(0);
        const Contour *c1 = sketch.contour(1);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                if (c0->curve(i) == c1->curve(j)) {
                    EXPECT_TRUE(c0->orientation(i) != c1->orientation(j));
                }
            }
        }
    }

    // Test Boundary
    {
        const Contour *boundary = sketch.boundary();

        EXPECT_EQ(boundary->size(), 4);
        EXPECT_TRUE(&l0 == boundary->curve(0) || &l0 == boundary->curve(1) || &l0 == boundary->curve(2) || &l0 == boundary->curve(3));
        EXPECT_TRUE(&l1 == boundary->curve(0) || &l1 == boundary->curve(1) || &l1 == boundary->curve(2) || &l1 == boundary->curve(3));
        EXPECT_TRUE(&l3 == boundary->curve(0) || &l3 == boundary->curve(1) || &l3 == boundary->curve(2) || &l3 == boundary->curve(3));
        EXPECT_TRUE(&l4 == boundary->curve(0) || &l4 == boundary->curve(1) || &l4 == boundary->curve(2) || &l4 == boundary->curve(3));
    }
}

TEST(Star, find_closed_contour_3) {
    /*
        Test contour with only two boundary curves
    */
    Sketch sketch;

    Vertex &vc = sketch.new_element<Vertex>(0.0, 0.0);
    Vertex &v0 = sketch.new_element<Vertex>(1.0, 0.0);
    Vertex &v1 = sketch.new_element<Vertex>(-1.0, 0.0);

    CircularArc &arc = sketch.new_element<CircularArc>(v0, v1, vc, 1.0);
    LineSegment &line = sketch.new_element<LineSegment>(v1, v0);

    // Manual contour construction
    {
        vector<Star> stars;
        stars.push_back(Star{&v0, &sketch});
        stars.push_back(Star{&v1, &sketch});

        star_angle_sum_equals_2pi(stars);
    }

    // Sketch internal contour parsing
    {
        sketch.solve();
        sketch.build();

        EXPECT_TRUE(sketch.size_contours() == 1);
        EXPECT_TRUE(*sketch.contour(0) == *sketch.boundary());

        const Contour *contour = sketch.contour(0);

        EXPECT_TRUE(contour->size() == 2);
        EXPECT_TRUE(&arc == contour->curve(0) || &arc == contour->curve(1));
        EXPECT_TRUE(&line == contour->curve(0) || &line == contour->curve(1));
    }
}

TEST(Star, with_construction_lines) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v2 = s.new_element<Vertex>(1.0, 1.0);
    Vertex &v3 = s.new_element<Vertex>(0.0, 1.0);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
    LineSegment &l3 = s.new_element<LineSegment>(v3, v0);
    LineSegment &l4 = s.new_element<LineSegment>(v0, v2);
    l4.ForConstruction = true;

    s.solve();
    s.build();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "Star__with_construction_lines");

    EXPECT_TRUE(s.size_contours() == 1);
    EXPECT_TRUE(s.contour(0)->size() == 4);
    EXPECT_TRUE(s.boundary()->size() == 4);
    EXPECT_TRUE(*s.contour(0) == *s.boundary());
}