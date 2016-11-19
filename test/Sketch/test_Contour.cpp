#include "test_Sketch.hpp"

TEST(Contour, Triangle_CCW) {
    Sketch s;

    auto v0 = s.new_element_SHARED_PTR<Vertex>(1.0,0.0);
    auto v1 = s.new_element_SHARED_PTR<Vertex>(1.0,0.0);
    auto v2 = s.new_element_SHARED_PTR<Vertex>(0.0,1.0);

    auto l0 = s.new_element_SHARED_PTR<LineSegment>(v0,v1);
    auto l1 = s.new_element_SHARED_PTR<LineSegment>(v1,v2);
    auto l2 = s.new_element_SHARED_PTR<LineSegment>(v2,v0);

    std::vector<std::shared_ptr<Curve>> vc{l0,l1,l2};

    Contour cont{vc};

    EXPECT_TRUE(cont.curve(0) == l0);
    EXPECT_TRUE(cont.curve(1) == l1);
    EXPECT_TRUE(cont.curve(2) == l2);

    EXPECT_TRUE(cont.vertex(0) == v0);
    EXPECT_TRUE(cont.vertex(1) == v1);
    EXPECT_TRUE(cont.vertex(2) == v2);

    EXPECT_TRUE(cont == cont);
}

TEST(Contour, Triangle_CW) {
    Sketch s;

    auto v0 = s.new_element_SHARED_PTR<Vertex>(1.0,0.0);
    auto v1 = s.new_element_SHARED_PTR<Vertex>(1.0,0.0);
    auto v2 = s.new_element_SHARED_PTR<Vertex>(0.0,1.0);

    auto l0 = s.new_element_SHARED_PTR<LineSegment>(v1,v0);
    auto l1 = s.new_element_SHARED_PTR<LineSegment>(v2,v1);
    auto l2 = s.new_element_SHARED_PTR<LineSegment>(v0,v2);

    std::vector<std::shared_ptr<Curve>> vc{l0,l1,l2};

    Contour cont{vc};

    EXPECT_TRUE(cont.curve(0) == l0);
    EXPECT_TRUE(cont.curve(1) == l2);
    EXPECT_TRUE(cont.curve(2) == l1);

    EXPECT_TRUE(cont.vertex(0) == v1);
    EXPECT_TRUE(cont.vertex(1) == v0);
    EXPECT_TRUE(cont.vertex(2) == v2);

    EXPECT_TRUE(cont == cont);
}

TEST(Contour, Triangle) {
    Sketch s;

    auto v0 = s.new_element_SHARED_PTR<Vertex>(1.0,0.0);
    auto v1 = s.new_element_SHARED_PTR<Vertex>(1.0,0.0);
    auto v2 = s.new_element_SHARED_PTR<Vertex>(0.0,1.0);

    auto l0 = s.new_element_SHARED_PTR<LineSegment>(v1,v0);
    auto l1 = s.new_element_SHARED_PTR<LineSegment>(v1,v2);
    auto l2 = s.new_element_SHARED_PTR<LineSegment>(v2,v0);

    std::vector<std::shared_ptr<Curve>> vc{l2,l0,l1};

    Contour cont{vc};

    EXPECT_TRUE(cont.curve(0) == l2);
    EXPECT_TRUE(cont.curve(1) == l0);
    EXPECT_TRUE(cont.curve(2) == l1);

    EXPECT_TRUE(cont.vertex(0) == v2);
    EXPECT_TRUE(cont.vertex(1) == v0);
    EXPECT_TRUE(cont.vertex(2) == v1);

    EXPECT_TRUE(cont == cont);
}

TEST(Contour, Nonclosed_Failure) {
    Sketch s;

    auto v0 = s.new_element_SHARED_PTR<Vertex>(0.0,0.0);
    auto v1 = s.new_element_SHARED_PTR<Vertex>(1.0,1.0);

    auto l0 = s.new_element_SHARED_PTR<LineSegment>(v0,v1);

    std::vector<std::shared_ptr<Curve>> c{l0};

    EXPECT_ANY_THROW(Contour cont{c}); // Construction should fail since the contour is not closed
}

TEST(Contour, Disjoint_Failure) {
    Sketch s;

    auto v0 = s.new_element_SHARED_PTR<Vertex>(0.0,0.0);
    auto v1 = s.new_element_SHARED_PTR<Vertex>(1.0,1.0);
    auto v2 = s.new_element_SHARED_PTR<Vertex>(0.0,1.0);
    auto v3 = s.new_element_SHARED_PTR<Vertex>(2.0,2.0);
    auto v4 = s.new_element_SHARED_PTR<Vertex>(3.0,3.0);
    auto v5 = s.new_element_SHARED_PTR<Vertex>(2.0,3.0);

    auto l0 = s.new_element_SHARED_PTR<LineSegment>(v0,v1);
    auto l1 = s.new_element_SHARED_PTR<LineSegment>(v1,v2);
    auto l2 = s.new_element_SHARED_PTR<LineSegment>(v2,v0);
    auto l3 = s.new_element_SHARED_PTR<LineSegment>(v3,v4);
    auto l4 = s.new_element_SHARED_PTR<LineSegment>(v4,v5);
    auto l5 = s.new_element_SHARED_PTR<LineSegment>(v5,v3);

    std::vector<std::shared_ptr<Curve>> c{l0,l1,l2,l3,l4,l5};

    EXPECT_ANY_THROW(Contour cont{c}); // Construction should fail since contour is not simple
}

TEST(Contour, Implicit_Self_Intersection_Failure) {
    Sketch s;

    auto v0 = s.new_element_SHARED_PTR<Vertex>(1.0,1.0);
    auto v1 = s.new_element_SHARED_PTR<Vertex>(-1.0,-1.0);
    auto v2 = s.new_element_SHARED_PTR<Vertex>(-1.0,1.0);
    auto v3 = s.new_element_SHARED_PTR<Vertex>(1.0,-1.0);

    auto l0 = s.new_element_SHARED_PTR<LineSegment>(v0,v1);
    auto l1 = s.new_element_SHARED_PTR<LineSegment>(v1,v2);
    auto l2 = s.new_element_SHARED_PTR<LineSegment>(v2,v3);
    auto l3 = s.new_element_SHARED_PTR<LineSegment>(v3,v0);

    std::vector<std::shared_ptr<Curve>> c{l0,l1,l2,l3};

    EXPECT_ANY_THROW(Contour cont{c}); // Construction should fail since contour is self intersecting

    /*
        TODO: Use nurbs representation to approximate potential intersection point.
        TODO: Use newton's method to solve intersection problem with the generated initial guesses.
    */
}

TEST(Contour, Explicit_Self_Intersection_Failure) {
    Sketch s;

    auto v0 = s.new_element_SHARED_PTR<Vertex>(1.0,1.0);
    auto v1 = s.new_element_SHARED_PTR<Vertex>(-1.0,-1.0);
    auto v2 = s.new_element_SHARED_PTR<Vertex>(-1.0,1.0);
    auto v3 = s.new_element_SHARED_PTR<Vertex>(1.0,-1.0);
    auto vc = s.new_element_SHARED_PTR<Vertex>(0.0,0.0);

    auto l0 = s.new_element_SHARED_PTR<LineSegment>(v0,vc);
    auto l1 = s.new_element_SHARED_PTR<LineSegment>(vc,v1);
    auto l2 = s.new_element_SHARED_PTR<LineSegment>(v1,v2);
    auto l3 = s.new_element_SHARED_PTR<LineSegment>(v2,vc);
    auto l4 = s.new_element_SHARED_PTR<LineSegment>(vc,v3);
    auto l5 = s.new_element_SHARED_PTR<LineSegment>(v3,v0);

    std::vector<std::shared_ptr<Curve>> c{l0,l1,l2,l3,l4,l5};

    EXPECT_ANY_THROW(Contour cont{c}); // Construction should fail since contour is self intersecting
}