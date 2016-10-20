#include "test_Sketch.hpp"

TEST(CONTOUR, CCW_TRIANGLE) {
    Vertex v0{1.0, 0.0};
    Vertex v1{1.0, 1.0};
    Vertex v2{0.0, 1.0};

    LineSegment l0{v0, v1};
    LineSegment l1{v1, v2};
    LineSegment l2{v2, v0};

    std::vector<const Curve *> vc(3);
    vc[0] = &l0;
    vc[1] = &l1;
    vc[2] = &l2;

    Contour cont{vc};

    EXPECT_TRUE(cont.curve(0) == &l0);
    EXPECT_TRUE(cont.curve(1) == &l1);
    EXPECT_TRUE(cont.curve(2) == &l2);

    EXPECT_TRUE(cont.vertex(0) == &v0);
    EXPECT_TRUE(cont.vertex(1) == &v1);
    EXPECT_TRUE(cont.vertex(2) == &v2);

    EXPECT_TRUE(cont == cont);
}

TEST(CONTOUR, CW_TRIANGLE) {
    Vertex v0{1.0, 0.0};
    Vertex v1{1.0, 1.0};
    Vertex v2{0.0, 1.0};

    LineSegment l0{v1, v0};
    LineSegment l1{v2, v1};
    LineSegment l2{v0, v2};

    std::vector<const Curve *> vc(3);
    vc[0] = &l0;
    vc[1] = &l1;
    vc[2] = &l2;

    Contour cont{vc};

    EXPECT_TRUE(cont.curve(0) == &l0);
    EXPECT_TRUE(cont.curve(1) == &l2);
    EXPECT_TRUE(cont.curve(2) == &l1);

    EXPECT_TRUE(cont.vertex(0) == &v1);
    EXPECT_TRUE(cont.vertex(1) == &v0);
    EXPECT_TRUE(cont.vertex(2) == &v2);

    EXPECT_TRUE(cont == cont);
}

TEST(CONTOUR, RANDOM_TRIANGLE) {
    Vertex v0{1.0, 0.0};
    Vertex v1{1.0, 1.0};
    Vertex v2{0.0, 1.0};

    LineSegment l0{v1, v0};
    LineSegment l1{v1, v2};
    LineSegment l2{v2, v0};

    std::vector<const Curve *> vc(3);
    vc[0] = &l2;
    vc[1] = &l0;
    vc[2] = &l1;

    Contour cont{vc};

    EXPECT_TRUE(cont.curve(0) == &l2);
    EXPECT_TRUE(cont.curve(1) == &l0);
    EXPECT_TRUE(cont.curve(2) == &l1);

    EXPECT_TRUE(cont.vertex(0) == &v2);
    EXPECT_TRUE(cont.vertex(1) == &v0);
    EXPECT_TRUE(cont.vertex(2) == &v1);

    EXPECT_TRUE(cont == cont);
}

TEST(CONTOUR, NONCLOSED_FAILURE) {
    Vertex v0{0.0, 0.0};
    Vertex v1{1.0, 1.0};
    LineSegment l0{v0, v1};
    std::vector<const Curve *> c(1);
    c[0] = &l0;

    EXPECT_ANY_THROW(Contour cont{c}); // Construction should fail since the contour is not closed
}

TEST(CONTOUR, DISJOINT_FAILURE) {
    Vertex v0{0.0, 0.0};
    Vertex v1{1.0, 1.0};
    Vertex v2{0.0, 1.0};
    LineSegment l0{v0, v1};
    LineSegment l1{v1, v2};
    LineSegment l2{v2, v0};

    Vertex v3{2.0, 2.0};
    Vertex v4{3.0, 3.0};
    Vertex v5{2.0, 3.0};
    LineSegment l3{v3, v4};
    LineSegment l4{v4, v5};
    LineSegment l5{v5, v3};

    std::vector<const Curve *> c(6);
    c[0] = &l0;
    c[1] = &l1;
    c[2] = &l2;
    c[3] = &l3;
    c[4] = &l4;
    c[5] = &l5;

    EXPECT_ANY_THROW(Contour cont{c}); // Construction should fail since contour is not simple
}

TEST(CONTOUR, IMPLICIT_SELF_INTERSECTION_FAILURE) {
    Vertex v0{1.0, 1.0};
    Vertex v1{-1.0, -1.0};
    Vertex v2{-1.0, 1.0};
    Vertex v3{1.0, -1.0};
    LineSegment l0{v0, v1};
    LineSegment l1{v1, v2};
    LineSegment l2{v2, v3};
    LineSegment l3{v3, v0};

    std::vector<const Curve *> c(4);
    c[0] = &l0;
    c[1] = &l1;
    c[2] = &l2;
    c[3] = &l3;

    EXPECT_ANY_THROW(Contour cont{c}); // Construction should fail since contour is self intersecting

    /*
        #TODO: Use nurbs representation to approximate potential intersection point.
        #TODO: Use newton's method to solve intersection problem with the generated initial guesses.
    */
}

TEST(CONTOUR, EXPLICIT_SELF_INTERSECTION_FAILURE) {
    Vertex v0{1.0, 1.0};
    Vertex v1{-1.0, -1.0};
    Vertex v2{-1.0, 1.0};
    Vertex v3{1.0, -1.0};
    Vertex vc{0.0, 0.0};

    LineSegment l0{v0, vc};
    LineSegment l1{vc, v1};
    LineSegment l2{v1, v2};
    LineSegment l3{v2, vc};
    LineSegment l4{vc, v3};
    LineSegment l5{v3, v0};

    std::vector<const Curve *> c(6);
    c[0] = &l0;
    c[1] = &l1;
    c[2] = &l2;
    c[3] = &l3;
    c[4] = &l4;
    c[5] = &l5;

    EXPECT_ANY_THROW(Contour cont{c}); // Construction should fail since contour is self intersecting
}