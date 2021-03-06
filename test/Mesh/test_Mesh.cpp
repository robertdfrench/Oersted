#include "test_Mesh.hpp"

TEST(Mesh, create_triangle_domain) {
    std::string test_name = "triangle_domain";

    Sketch s;

    auto v0 = s.new_element<Vertex>(1.0, 0.0);
    auto v1 = s.new_element<Vertex>(0.0, sqrt(3.0));
    auto v2 = s.new_element<Vertex>(-1.0, 0.0);

    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v1, v2);
    auto l2 = s.new_element<LineSegment>(v2, v0);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);
    bool result = s.build();
    ASSERT_TRUE(result);

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, test_name + "_mesh");

    std::vector<size_t> vmap = map_verticies_to_points({v0, v1, v2}, m);

    // Test number of vertices, edges, triangles
    EXPECT_TRUE(m.size_points() == 3);
    EXPECT_TRUE(m.size_edges() == 3);
    EXPECT_TRUE(m.size_triangles() == 1);

    EXPECT_TRUE(m.num_points() == 3);
    EXPECT_TRUE(m.num_edges() == 3);
    EXPECT_TRUE(m.num_triangles() == 1);

    EXPECT_TRUE(v0->x() == m.point(vmap[0]).X);
    EXPECT_TRUE(v0->y() == m.point(vmap[0]).Y);

    EXPECT_TRUE(v1->x() == m.point(vmap[1]).X);
    EXPECT_TRUE(v1->y() == m.point(vmap[1]).Y);

    EXPECT_TRUE(v2->x() == m.point(vmap[2]).X);
    EXPECT_TRUE(v2->y() == m.point(vmap[2]).Y);

    // Test validity, optimality
    edges_are_valid(m);
    edges_are_optimal(m);

    // Test edge and node connections
    const Edge e0 = m.edge(vmap[0]);
    const Edge e1 = m.edge(vmap[1]);
    const Edge e2 = m.edge(vmap[2]);

    EXPECT_TRUE(v0->x() == m.point(e0.node()).X);
    EXPECT_TRUE(v0->y() == m.point(e0.node()).Y);
    EXPECT_TRUE(e1.self() == e0.next());
    EXPECT_TRUE(e2.self() == e0.prev());

    EXPECT_TRUE(v1->x() == m.point(e1.node()).X);
    EXPECT_TRUE(v1->y() == m.point(e1.node()).Y);
    EXPECT_TRUE(e2.self() == e1.next());
    EXPECT_TRUE(e0.self() == e1.prev());

    EXPECT_TRUE(v2->x() == m.point(e2.node()).X);
    EXPECT_TRUE(v2->y() == m.point(e2.node()).Y);
    EXPECT_TRUE(e0.self() == e2.next());
    EXPECT_TRUE(e1.self() == e2.prev());

    // Test triangles
    const Edge e = m.triangle(0);
    Point cc = m.circumcenter(e.self());
    EXPECT_NEAR(0.0, cc.X, TOL);
    EXPECT_NEAR(sqrt(3.0) / 3.0, cc.Y, TOL);

    double cr = m.circumradius(e.self());
    EXPECT_NEAR(2.0 * sqrt(3.0) / 3.0, cr, TOL);

    // Forced Refinement
    forced_refinement(m, test_name + "_mesh_refine_loop", 7);

    // Test refinement algorithm
    m = Mesh(s);
    m.create();
    m.MaximumElementSize = 0.1;
    m.MinimumElementSize = 0.01;
    m.MinimumElementQuality = 0.1;
    m.refine();
    m.save_as(SAVE_DIR, test_name + "_mesh_refine_algorithm");
}

TEST(Mesh, create_square_domain) {
    std::string test_name = "square_domain";

    Sketch s;

    auto v0 = s.new_element<Vertex>(0.0, 0.0);
    auto v1 = s.new_element<Vertex>(1.0, 0.0);
    auto v2 = s.new_element<Vertex>(1.0, 1.0);
    auto v3 = s.new_element<Vertex>(0.0, 1.0);

    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v1, v2);
    auto l2 = s.new_element<LineSegment>(v2, v3);
    auto l3 = s.new_element<LineSegment>(v3, v0);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);
    bool result = s.build();
    ASSERT_TRUE(result);

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, test_name + "_mesh");

    // Test number of verticies, edges, triangles
    std::vector<size_t> vmap = map_verticies_to_points({v0, v1, v2, v3}, m);

    EXPECT_TRUE(m.size_points() == 4);
    EXPECT_TRUE(m.size_edges() == 6);
    EXPECT_TRUE(m.size_triangles() == 2);

    EXPECT_TRUE(m.num_points() == 4);
    EXPECT_TRUE(m.num_edges() == 5);
    EXPECT_TRUE(m.num_triangles() == 2);

    EXPECT_TRUE(v0->x() == m.point(vmap[0]).X);
    EXPECT_TRUE(v0->y() == m.point(vmap[0]).Y);

    EXPECT_TRUE(v1->x() == m.point(vmap[1]).X);
    EXPECT_TRUE(v1->y() == m.point(vmap[1]).Y);

    EXPECT_TRUE(v2->x() == m.point(vmap[2]).X);
    EXPECT_TRUE(v2->y() == m.point(vmap[2]).Y);

    EXPECT_TRUE(v3->x() == m.point(vmap[3]).X);
    EXPECT_TRUE(v3->y() == m.point(vmap[3]).Y);

    // Test validity, optimality
    edges_are_valid(m);
    edges_are_optimal(m);

    // Test edge and node connections
    EXPECT_TRUE(m.point((size_t) 0) == m.point(m.edge(0)));
    EXPECT_TRUE(m.point((size_t) 1) == m.point(m.edge(1)));
    EXPECT_TRUE(m.point((size_t) 2) == m.point(m.edge(2)));
    EXPECT_TRUE(m.point((size_t) 3) == m.point(m.edge(3)));
    EXPECT_TRUE(m.point((size_t) 1) == m.point(m.edge(4)));
    EXPECT_TRUE(m.point((size_t) 3) == m.point(m.edge(5)));

    for (size_t i = 0; i < 5; i++) {
        const Edge e = m.edge(i);
        EXPECT_TRUE(e.self() == m.edge(e.next()).prev());
        EXPECT_TRUE(e.self() == m.edge(e.prev()).next());
    }

    // Test triangles
    for (size_t i = 0; i < m.size_triangles(); ++i) {
        Point cc = m.circumcenter(m.triangle(0).self());
        EXPECT_NEAR(0.5, cc.X, TOL);
        EXPECT_NEAR(0.5, cc.Y, TOL);
    }

    // Forced Refinement
    forced_refinement(m, test_name + "_mesh_refine_loop", 7);

    // Test refinement algorithm
    m = Mesh(s);
    m.create();
    m.MaximumElementSize = 0.1;
    m.MinimumElementSize = 0.01;
    m.MinimumElementQuality = 0.1;
    m.refine();
    m.save_as(SAVE_DIR, test_name + "_mesh_refine_algorithm");
}

TEST(Mesh, create_narrow_diamond_domain) {
    std::string test_name = "narrow_diamond_domain";

    Sketch s;

    auto v0 = s.new_element<Vertex>(1.0, 0.0);
    auto v1 = s.new_element<Vertex>(0.0, 2.0);
    auto v2 = s.new_element<Vertex>(-1.0, 0.0);
    auto v3 = s.new_element<Vertex>(0.0, -2.0);

    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v1, v2);
    auto l2 = s.new_element<LineSegment>(v2, v3);
    auto l3 = s.new_element<LineSegment>(v3, v0);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);
    bool result = s.build();
    ASSERT_TRUE(result);

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, test_name + "_mesh");

    // Test number of verticies, edges, triangles
    EXPECT_TRUE(m.size_points() == 4);
    EXPECT_TRUE(m.size_edges() == 6);
    EXPECT_TRUE(m.size_triangles() == 2);

    EXPECT_TRUE(m.num_points() == 4);
    EXPECT_TRUE(m.num_edges() == 5);
    EXPECT_TRUE(m.num_triangles() == 2);

    // Test validity, optimality
    edges_are_valid(m);
    edges_are_optimal(m);

    // Test for proper edge swaps
    std::vector<size_t> vmap = map_verticies_to_points({v0, v1, v2, v3}, m);

    for (size_t i = 5; i < 6; i++) {
        const Edge e = m.edge(5);

        if (m.point(e.node()) == m.point(vmap[0])) {
            EXPECT_TRUE(m.point(vmap[2]) == m.point(m.edge(e.next()).node()));
            EXPECT_TRUE(m.point(vmap[2]) == m.point(m.edge(e.twin()).node()));
        } else if (m.point(e.node()) == m.point(vmap[2])) {
            EXPECT_TRUE(m.point(vmap[0]) == m.point(m.edge(e.next()).node()));
            EXPECT_TRUE(m.point(vmap[0]) == m.point(m.edge(e.twin()).node()));
        }
    }

    // Test triangle circumcenters
    Point cc0 = m.circumcenter(m.triangle(0).self());
    Point cc1 = m.circumcenter(m.triangle(1).self());

    EXPECT_NEAR(0.0, cc0.X, TOL);
    EXPECT_NEAR(0.75, std::abs(cc0.Y), TOL);

    EXPECT_NEAR(0.0, cc1.X, TOL);
    EXPECT_NEAR(0.75, std::abs(cc1.Y), TOL);

    EXPECT_NEAR(-cc0.Y, cc1.Y, TOL);
    // Forced Refinement
    forced_refinement(m, test_name + "_mesh_refine_loop", 7);

    // Test refinement algorithm
    m = Mesh(s);
    m.create();
    m.MaximumElementSize = 0.1;
    m.MinimumElementSize = 0.01;
    m.MinimumElementQuality = 0.1;
    m.refine();
    m.save_as(SAVE_DIR, test_name + "_mesh_refine_algorithm");
}

TEST(Mesh, create_narrow_rectangle_domain) {
    std::string test_name = "narrow_rectangle_domain";

    Sketch s;

    auto v0 = s.new_element<Vertex>(5.0, 0.0);
    auto v1 = s.new_element<Vertex>(5.0, 1.0);
    auto v2 = s.new_element<Vertex>(0.0, 2.0);
    auto v3 = s.new_element<Vertex>(-5.0, 1.0);
    auto v4 = s.new_element<Vertex>(-5.0, 0.0);

    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v1, v2);
    auto l2 = s.new_element<LineSegment>(v2, v3);
    auto l3 = s.new_element<LineSegment>(v3, v4);
    auto l4 = s.new_element<LineSegment>(v4, v0);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    bool result = s.build();
    ASSERT_TRUE(result);

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, test_name + "_mesh");

    // Test number of verticies, edges, triangles
    EXPECT_TRUE(m.size_points() == 6);
    EXPECT_TRUE(m.size_edges() == 12);
    EXPECT_TRUE(m.size_triangles() == 4);

    EXPECT_TRUE(m.num_points() == 6);
    EXPECT_TRUE(m.num_edges() == 9);
    EXPECT_TRUE(m.num_triangles() == 4);

    // Test validity, optimality
    edges_are_valid(m);
    edges_are_optimal(m);

    // Test edge splits
    std::vector<size_t> vmap = map_verticies_to_points({v0, v1, v2, v3, v4}, m);

    Point const &v5 = m.point(5);
    EXPECT_NEAR(0.0, v5.X, DBL_EPSILON);
    EXPECT_NEAR(0.0, v5.Y, DBL_EPSILON);

    for (size_t i = 0; i < 12; ++i) {
        Edge const e = m.edge(i);

        if (m.point(e.node()) == m.point(vmap[2])) {
            if (e.twin() != e.self()) {
                EXPECT_TRUE(m.point(m.edge(e.twin()).node()) == v5);
            }
        }
    }

    // Forced Refinement
    forced_refinement(m, test_name + "_mesh_refine_loop", 7);

    // Test refinement algorithm
    m = Mesh(s);
    m.create();
    m.MaximumElementSize = 0.1;
    m.MinimumElementSize = 0.01;
    m.MinimumElementQuality = 0.1;
    m.refine();
    m.save_as(SAVE_DIR, test_name + "_mesh_refine_algorithm");
}

TEST(Mesh, create_half_circle_domain) {
    std::string test_name = "half_circle_domain";

    Sketch s;

    auto vc = s.new_element<Vertex>(0.0, 0.0);
    auto v0 = s.new_element<Vertex>(1.0, 0.0);
    auto v1 = s.new_element<Vertex>(-1.0, 0.0);

    auto arc = s.new_element<CircularArc>(v0, v1, vc, 1.0);
    auto line = s.new_element<LineSegment>(v1, v0);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);
    bool result = s.build();
    ASSERT_TRUE(result);

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, test_name + "_mesh");

    // Test number of verticies, edges, triangles
    EXPECT_TRUE(m.size_points() == 6);
    EXPECT_TRUE(m.size_edges() == 12);
    EXPECT_TRUE(m.size_triangles() == 4);

    EXPECT_TRUE(m.num_points() == 6);
    EXPECT_TRUE(m.num_edges() == 9);
    EXPECT_TRUE(m.num_triangles() == 4);

    // Test validity, optimality
    edges_are_valid(m);
    edges_are_optimal(m);

    // Forced Refinement
    forced_refinement(m, test_name + "_mesh_refine_loop", 7);

    // Test refinement algorithm
    m = Mesh(s);
    m.create();
    m.MaximumElementSize = 0.1;
    m.MinimumElementSize = 0.01;
    m.MinimumElementQuality = 0.1;
    m.refine();
    m.save_as(SAVE_DIR, test_name + "_mesh_refine_algorithm");
}

TEST(Mesh, create_horseshoe_domain) {
    std::string test_name{"horseshoe_domain"};

    Sketch s;

    auto vc = s.new_element<Vertex>(0.0, 0.0);
    auto v0 = s.new_element<Vertex>(2.0, 0.0);
    auto v1 = s.new_element<Vertex>(-2.0, 0.0);
    auto v2 = s.new_element<Vertex>(-1.9, 0.0);
    auto v3 = s.new_element<Vertex>(1.9, 0.0);

    auto arc0 = s.new_element<CircularArc>(v0, v1, vc, 2.0);
    auto line0 = s.new_element<LineSegment>(v1, v2);
    auto arc1 = s.new_element<CircularArc>(v3, v2, vc, 1.9);
    auto line1 = s.new_element<LineSegment>(v3, v0);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);
    bool result = s.build();
    ASSERT_TRUE(result);

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, test_name + "_mesh");

    { // Test triangles, possibly redundant
        for (size_t i = 0; i < m.size_edges(); ++i) {
            EXPECT_TRUE(m.edge(m.edge(i).next()).prev() == m.edge(i).self());
            EXPECT_TRUE(m.edge(m.edge(m.edge(m.edge(i).next()).next()).prev()).prev() == m.edge(i).self());
            EXPECT_TRUE(m.edge(m.edge(m.edge(i).next()).next()).next() == m.edge(i).self());

            EXPECT_TRUE(m.edge(m.edge(i).prev()).next() == m.edge(i).self());
            EXPECT_TRUE(m.edge(m.edge(m.edge(m.edge(i).prev()).prev()).next()).next() == m.edge(i).self());
            EXPECT_TRUE(m.edge(m.edge(m.edge(i).prev()).prev()).prev() == m.edge(i).self());

            if ((m.edge(i).twin() != m.edge(i).self())) {
                EXPECT_TRUE(m.edge(m.edge(i).twin()).node() == m.edge(m.edge(i).next()).node());
            }
        }
    }

    // Test validity, optimality
    edges_are_valid(m);
    edges_are_optimal(m);

    // Forced refinement
    forced_refinement(m, test_name + "_mesh_refine_loop", 7);


    // Test refinement algorithm
    m = Mesh(s);
    m.create();
    m.MaximumElementSize = 0.1;
    m.MinimumElementSize = 0.01;
    m.MinimumElementQuality = 0.1;
    m.refine();
    m.save_as(SAVE_DIR, test_name + "_mesh_refine_algorithm");
}

TEST(Mesh, create_I_shaped_domain) {
    std::string test_name{"i_domain"};

    Sketch s;

    auto v0 = s.new_element<Vertex>(0.0, 0.0);
    auto v1 = s.new_element<Vertex>(3.0, 0.0);
    auto v2 = s.new_element<Vertex>(3.0, 1.0);
    auto v3 = s.new_element<Vertex>(2.0, 1.0);
    auto v4 = s.new_element<Vertex>(2.0, 2.0);
    auto v5 = s.new_element<Vertex>(3.0, 2.0);
    auto v6 = s.new_element<Vertex>(3.0, 3.0);
    auto v7 = s.new_element<Vertex>(0.0, 3.0);
    auto v8 = s.new_element<Vertex>(0.0, 2.0);
    auto v9 = s.new_element<Vertex>(1.0, 2.0);
    auto v10 = s.new_element<Vertex>(1.0, 1.0);
    auto v11 = s.new_element<Vertex>(0.0, 1.0);

    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v1, v2);
    auto l2 = s.new_element<LineSegment>(v2, v3);
    auto l3 = s.new_element<LineSegment>(v3, v4);
    auto l4 = s.new_element<LineSegment>(v4, v5);
    auto l5 = s.new_element<LineSegment>(v5, v6);
    auto l6 = s.new_element<LineSegment>(v6, v7);
    auto l7 = s.new_element<LineSegment>(v7, v8);
    auto l8 = s.new_element<LineSegment>(v8, v9);
    auto l9 = s.new_element<LineSegment>(v9, v10);
    auto l10 = s.new_element<LineSegment>(v10, v11);
    auto l11 = s.new_element<LineSegment>(v11, v0);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, test_name);

    bool result = s.build();
    ASSERT_TRUE(result);

    EXPECT_TRUE(s.size_contours() == 1);

    Mesh m{s};

    m.create();

    m.save_as(SAVE_DIR, test_name + "_mesh");

    // Test validity, optimality
    edges_are_valid(m);
    edges_are_optimal(m);

    // Forced Refinement
    forced_refinement(m, test_name + "_mesh_refine_loop", 7);

    // Test refinement algorithm
    m = Mesh(s);
    m.create();
    m.MaximumElementSize = 0.1;
    m.MinimumElementSize = 0.01;
    m.MinimumElementQuality = 0.1;
    m.refine();
    m.save_as(SAVE_DIR, test_name + "_mesh_refine_algorithm");
}

TEST(Mesh, create_corner_square_domain) {
    std::string test_name{"corner_square_domain"};

    Sketch s;

    auto v0 = s.new_element<Vertex>(0.0, 0.0);
    auto v1 = s.new_element<Vertex>(0.5, 0.0);
    auto v2 = s.new_element<Vertex>(1.0, 0.0);
    auto v3 = s.new_element<Vertex>(1.0, 0.5);
    auto v4 = s.new_element<Vertex>(1.0, 1.0);
    auto v5 = s.new_element<Vertex>(0.0, 1.0);
    auto v6 = s.new_element<Vertex>(0.5, 0.5);

    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v1, v6);
    auto l2 = s.new_element<LineSegment>(v6, v3);
    auto l3 = s.new_element<LineSegment>(v3, v4);
    auto l4 = s.new_element<LineSegment>(v4, v5);
    auto l5 = s.new_element<LineSegment>(v5, v0);

    auto l6 = s.new_element<LineSegment>(v1, v2);
    auto l7 = s.new_element<LineSegment>(v2, v3);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);

    bool result = s.build();
    ASSERT_TRUE(result);

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, test_name + "_mesh");

    // Test number of vertices, edges, triangles
    EXPECT_TRUE(m.num_points() == 9);
    EXPECT_TRUE(m.size_points() == 9);

    EXPECT_TRUE(m.num_edges() == 16);
    EXPECT_TRUE(m.size_edges() == 24);

    EXPECT_TRUE(m.num_triangles() == 8);
    EXPECT_TRUE(m.size_triangles() == 8);

    // Test validity, optimality
    edges_are_valid(m);
    edges_are_optimal(m);

    // Forced refinement
    forced_refinement(m, test_name + "_mesh_refine_loop", 7);

    // Test refinement algorithm
    m = Mesh(s);
    m.create();
    m.MaximumElementSize = 0.1;
    m.MinimumElementSize = 0.01;
    m.MinimumElementQuality = 0.1;
    m.refine();
    m.save_as(SAVE_DIR, test_name + "_mesh_refine_algorithm");
}

TEST(Mesh, create_square_in_square_domain) {
    std::string test_name{"square_in_square"};

    Sketch s;

    auto v0 = s.new_element<Vertex>(0.0, 0.0);
    auto v1 = s.new_element<Vertex>(2.0, 0.0);
    auto v2 = s.new_element<Vertex>(2.0, 2.0);
    auto v3 = s.new_element<Vertex>(0.0, 2.0);

    auto v4 = s.new_element<Vertex>(0.9, 0.9);
    auto v5 = s.new_element<Vertex>(1.9, 0.9);
    auto v6 = s.new_element<Vertex>(1.9, 1.9);
    auto v7 = s.new_element<Vertex>(0.9, 1.9);

    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v1, v2);
    auto l2 = s.new_element<LineSegment>(v2, v3);
    auto l3 = s.new_element<LineSegment>(v3, v0);

    auto l4 = s.new_element<LineSegment>(v4, v5);
    auto l5 = s.new_element<LineSegment>(v5, v6);
    auto l6 = s.new_element<LineSegment>(v6, v7);
    auto l7 = s.new_element<LineSegment>(v7, v4);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);
    bool result = s.build();
    ASSERT_TRUE(result);

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, test_name);

    Mesh m{s};
    m.MaximumElementSize = 0.1;
    m.MinimumElementSize = 0.01;
    m.MinimumElementQuality = 0.1;

    m.create();

    m.save_as(SAVE_DIR, test_name + "_mesh");

    edges_are_valid(m);
    edges_are_optimal(m);

    m.refine();

    m.save_as(SAVE_DIR, test_name + "_mesh_refine_algorithm");
}

/*
TEST(Edge, are_intersecting__nonparallel_true) {
    Point v0{0.0, 0.0};
    Point v1{1.0, 1.0};
    Point v2{0.0, 1.0};
    Point v3{1.0, 0.0};

    Edge dummy;
    Edge e1{v1, dummy, dummy, dummy};
    Edge e0{v0, e1, dummy, dummy};

    Edge e3{v3, dummy, dummy, dummy};
    Edge e2{v2, e3, dummy, dummy};

    EXPECT_TRUE(are_intersecting(&e0, &e2));
}

TEST(Edge, are_intersecting__nonparallel_false) {
    Point v0{1.1, 1.1};
    Point v1{2.0, 2.0};
    Point v2{0.0, 1.0};
    Point v3{1.0, 0.0};

    Edge dummy;
    Edge e1{v1, dummy, dummy, dummy};
    Edge e0{v0, e1, dummy, dummy};

    Edge e3{v3, dummy, dummy, dummy};
    Edge e2{v2, e3, dummy, dummy};

    EXPECT_FALSE(are_intersecting(&e0, &e2));
}

TEST(Edge, are_intersecting__nonparallel_false_shared_vertex) {
    Point v0{0.0, 0.0};
    Point v1{1.0, 1.0};
    Point v2{0.0, 2.0};

    Edge dummy;
    Edge e1{v1, dummy, dummy, dummy};
    Edge e0{v0, e1, dummy, dummy};

    Edge e3{v2, dummy, dummy, dummy};
    Edge e2{v1, e3, dummy, dummy};

    EXPECT_FALSE(are_intersecting(&e0, &e2));
}

TEST(Edge, are_intersecting__nonparallel_false_coincident_vertex) {
    Point v0{0.0, 0.0};
    Point v1{1.0, 1.0};
    Point v2{0.0, 1.0};
    Point v3{0.5, 0.5};

    Edge dummy;
    Edge e1{v1, dummy, dummy, dummy};
    Edge e0{v0, e1, dummy, dummy};

    Edge e3{v3, dummy, dummy, dummy};
    Edge e2{v2, e3, dummy, dummy};

    EXPECT_FALSE(are_intersecting(&e0, &e2));
}

TEST(Edge, are_intersecting__parallel_shaded_false) {
    Point v0{0.0, 0.0};
    Point v1{1.0, 1.0};
    Point v2{0.0, 0.5};
    Point v3{1.0, 1.5};

    Edge dummy;
    Edge e1{v1, dummy, dummy, dummy};
    Edge e0{v0, e1, dummy, dummy};

    Edge e3{v3, dummy, dummy, dummy};
    Edge e2{v2, e3, dummy, dummy};

    EXPECT_FALSE(are_intersecting(&e0, &e2));
}

TEST(Edge, are_intersecting__parallel_nonshaded_false) {
    Point v0{0.0, 0.0};
    Point v1{1.0, 1.0};
    Point v2{0.0, 1.5};
    Point v3{1.0, 2.5};

    Edge dummy;
    Edge e1{v1, dummy, dummy, dummy};
    Edge e0{v0, e1, dummy, dummy};

    Edge e3{v3, dummy, dummy, dummy};
    Edge e2{v2, e3, dummy, dummy};

    EXPECT_FALSE(are_intersecting(&e0, &e2));
}

TEST(Edge, are_intersecting__colinear_true) {
    Point v0{0.0, 0.0};
    Point v1{1.0, 1.0};
    Point v2{0.5, 0.5};
    Point v3{1.5, 1.5};

    Edge dummy;
    Edge e1{v1, dummy, dummy, dummy};
    Edge e0{v0, e1, dummy, dummy};

    Edge e3{v3, dummy, dummy, dummy};
    Edge e2{v2, e3, dummy, dummy};

    EXPECT_TRUE(are_intersecting(&e0, &e2));
}

TEST(Edge, are_intersecting__colinear_false) {
    Point v0{0.0, 0.0};
    Point v1{1.0, 1.0};
    Point v2{1.5, 1.5};
    Point v3{2.5, 2.5};

    Edge dummy;
    Edge e1{v1, dummy, dummy, dummy};
    Edge e0{v0, e1, dummy, dummy};

    Edge e3{v3, dummy, dummy, dummy};
    Edge e2{v2, e3, dummy, dummy};

    EXPECT_FALSE(are_intersecting(&e0, &e2));
}

TEST(Edge, are_intersecting__colinear_false_shared_vertex) {
    Point v0{0.0, 0.0};
    Point v1{1.0, 1.0};
    Point v2{2.0, 2.0};

    Edge dummy;
    Edge e1{v1, dummy, dummy, dummy};
    Edge e0{v0, e1, dummy, dummy};

    Edge e3{v2, dummy, dummy, dummy};
    Edge e2{v1, e3, dummy, dummy};

    EXPECT_FALSE(are_intersecting(&e0, &e2));
}
*/

TEST(Mesh, locate_triangle__triangular_domain) {
    Sketch s;

    auto v0 = s.new_element<Vertex>(0.0, 0.0);
    auto v1 = s.new_element<Vertex>(2.0, 0.0);
    auto v2 = s.new_element<Vertex>(1.0, sqrt(3.0));

    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v1, v2);
    auto l2 = s.new_element<LineSegment>(v2, v0);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);
    bool result = s.build();
    ASSERT_TRUE(result);

    /*
    std::vector<const Curve*> cc{ &l0,&l1,&l2 };
    Contour cont{ cc };

    std::vector<const Contour*> cv{ &cont };
    Mesh mesh{ cv };
    */
    Mesh mesh{s};
    mesh.create();

    EXPECT_TRUE(mesh.size_edges() == 3);
    EXPECT_TRUE(mesh.size_points() == 3);

    // Interior and Exterior
    Point vi{1.0, 1.0};
    Point ve0{1.0, -1.0};
    Point ve1{2.0, 2.0};
    Point ve2{0.0, 2.0};

    std::vector<size_t> vmap = map_verticies_to_points({v0, v1, v2}, mesh);

    for (size_t i = 0; i < mesh.size_edges(); ++i) {
        Edge e;
        Point vp;
        size_t loc = i;

        vp = vi;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Interior);

        vp = ve0;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[0])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[1])));

        vp = ve1;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[1])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[2])));

        vp = ve2;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[2])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[0])));
    }

    // Corner Cases
    ve0 = {1.0, -DBL_EPSILON};
    ve1 = {1.5 + DBL_EPSILON, sqrt(3.0) / 2.0 + DBL_EPSILON};
    ve2 = {0.5 - DBL_EPSILON, sqrt(3.0) / 2.0 + DBL_EPSILON};

    for (size_t i = 0; i != mesh.size_edges(); ++i) {
        Edge e;
        Point vp;
        size_t loc = i;

        vp = ve0;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[0])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[1])));

        vp = ve1;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[1])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[2])));

        vp = ve2;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[2])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[0])));
    }
}

TEST(Mesh, locate_triange__square_domain) {
    Sketch s;

    auto v0 = s.new_element<Vertex>(0.0, 0.0);
    auto v1 = s.new_element<Vertex>(1.0, 0.0);
    auto v2 = s.new_element<Vertex>(1.0, 1.0);
    auto v3 = s.new_element<Vertex>(0.0, 1.0);
    auto l0 = s.new_element<LineSegment>(v0, v1);
    auto l1 = s.new_element<LineSegment>(v1, v2);
    auto l2 = s.new_element<LineSegment>(v2, v3);
    auto l3 = s.new_element<LineSegment>(v3, v0);

    double res_norm = s.solve();
    EXPECT_LE(res_norm, FLT_EPSILON);
    bool build_result = s.build();
    ASSERT_TRUE(build_result);

    Mesh mesh{s};
    mesh.create();

    EXPECT_TRUE(mesh.size_edges() == 6);
    EXPECT_TRUE(mesh.size_points() == 4);

    // Interior Points
    Point vi0{0.25, 0.25};
    Point vi1{0.75, 0.75};
    for (size_t i = 0; i < mesh.size_edges(); ++i) {
        Point vp;
        size_t loc = i;

        vp = vi0;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Interior);

        vp = vi1;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Interior);
    }

    // Exterior Points
    std::vector<size_t> vmap = map_verticies_to_points({v0, v1, v2, v3}, mesh);

    Point ve0{0.5, -1.0};
    Point ve1{2.0, 0.5};
    Point ve2{0.5, 2.0};
    Point ve3{-1.0, 0.5};
    for (size_t i = 0; i < mesh.size_edges(); ++i) {
        Edge e;
        Point vp;
        size_t loc = i;

        vp = ve0;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[0])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[1])));

        vp = ve1;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[1])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[2])));

        vp = ve2;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[2])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[3])));

        vp = ve3;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[3])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[0])));
    }

    // Exterior Edge Points
    ve0 = {0.5, -DBL_EPSILON};
    ve1 = {1.0 + DBL_EPSILON, 0.5};
    ve2 = {0.5, 1.0 + DBL_EPSILON};
    ve3 = {0.0 - DBL_EPSILON, 0.5};
    LocateTriangleResult result;
    for (size_t i = 0; i < mesh.size_edges(); ++i) {
        Edge e;
        Point vp;
        size_t loc = i;

        vp = ve0;
        result = mesh.locate_triangle(vp, loc);
        EXPECT_TRUE(result == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[0])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[1])));

        vp = ve1;
        result = mesh.locate_triangle(vp, loc);
        EXPECT_TRUE(result == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[1])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[2])));

        vp = ve2;
        result = mesh.locate_triangle(vp, loc);
        EXPECT_TRUE(result == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[2])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[3])));

        vp = ve3;
        result = mesh.locate_triangle(vp, loc);
        EXPECT_TRUE(result == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e.node()) == mesh.point(vmap[3])) && (mesh.point(mesh.edge(e.next()).node()) == mesh.point(vmap[0])));
    }

    Point vie{0.5, 0.5};
    for (size_t i = 0; i < mesh.size_edges(); ++i) {
        EXPECT_TRUE(mesh.locate_triangle(vie, i) == LocateTriangleResult::Interior);
    }
}