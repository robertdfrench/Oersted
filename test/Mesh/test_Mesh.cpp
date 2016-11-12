#include "test_Mesh.hpp"

TEST(Mesh, create__triangle_domain) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(0.0, sqrt(3.0));
    Vertex &v2 = s.new_element<Vertex>(-1.0, 0.0);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l2 = s.new_element<LineSegment>(v2, v0);

    s.solve();
    s.build();

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, "triangle_domain");

    std::vector<size_t> vmap = map_verticies_to_points({v0, v1, v2}, m);

    { // Test number of vertices, edges, triangles
        EXPECT_TRUE(m.size_points() == 3);
        EXPECT_TRUE(m.size_edges() == 3);
        EXPECT_TRUE(m.size_triangles() == 1);

        EXPECT_TRUE(m.num_points() == 3);
        EXPECT_TRUE(m.num_edges() == 3);
        EXPECT_TRUE(m.num_triangles() == 1);

        EXPECT_TRUE(v0.x() == m.point(vmap[0]).X);
        EXPECT_TRUE(v0.y() == m.point(vmap[0]).Y);

        EXPECT_TRUE(v1.x() == m.point(vmap[1]).X);
        EXPECT_TRUE(v1.y() == m.point(vmap[1]).Y);

        EXPECT_TRUE(v2.x() == m.point(vmap[2]).X);
        EXPECT_TRUE(v2.y() == m.point(vmap[2]).Y);
    }

    { // Test validity, optimality
        edges_are_valid(m);
        edges_are_optimal(m);
    }


    { // Test edge and node connections
        const Edge *e0 = m.edge(vmap[0]);
        const Edge *e1 = m.edge(vmap[1]);
        const Edge *e2 = m.edge(vmap[2]);

        EXPECT_TRUE(v0.x() == m.point(e0->node()).X);
        EXPECT_TRUE(v0.y() == m.point(e0->node()).Y);
        EXPECT_TRUE(e1->self() == e0->next());
        EXPECT_TRUE(e2->self() == e0->prev());

        EXPECT_TRUE(v1.x() == m.point(e1->node()).X);
        EXPECT_TRUE(v1.y() == m.point(e1->node()).Y);
        EXPECT_TRUE(e2->self() == e1->next());
        EXPECT_TRUE(e0->self() == e1->prev());

        EXPECT_TRUE(v2.x() == m.point(e2->node()).X);
        EXPECT_TRUE(v2.y() == m.point(e2->node()).Y);
        EXPECT_TRUE(e0->self() == e2->next());
        EXPECT_TRUE(e1->self() == e2->prev());
    }

    { // Test triangles
        const Edge *e = m.triangle(0);
        Point cc = m.circumcenter(e->self());
        EXPECT_NEAR(0.0, cc.X, TOL);
        EXPECT_NEAR(sqrt(3.0) / 3.0, cc.Y, TOL);

        double cr = m.circumradius(e->self());
        EXPECT_NEAR(2.0 * sqrt(3.0) / 3.0, cr, TOL);
    }

    { // Forced Refinement
        forced_refinement(m, "triangle_domain_refine_loop", 7);
    }

    m.delete_me();
}

TEST(Mesh, create__square_domain) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v2 = s.new_element<Vertex>(1.0, 1.0);
    Vertex &v3 = s.new_element<Vertex>(0.0, 1.0);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
    LineSegment &l3 = s.new_element<LineSegment>(v3, v0);

    s.solve();
    s.build();

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, "square_domain");

    // Test number of verticies, edges, triangles
    {
        std::vector<size_t> vmap = map_verticies_to_points({v0, v1, v2, v3}, m);

        EXPECT_TRUE(m.size_points() == 4);
        EXPECT_TRUE(m.size_edges() == 6);
        EXPECT_TRUE(m.size_triangles() == 2);

        EXPECT_TRUE(m.num_points() == 4);
        EXPECT_TRUE(m.num_edges() == 5);
        EXPECT_TRUE(m.num_triangles() == 2);

        EXPECT_TRUE(v0.x() == m.point(vmap[0]).X);
        EXPECT_TRUE(v0.y() == m.point(vmap[0]).Y);

        EXPECT_TRUE(v1.x() == m.point(vmap[1]).X);
        EXPECT_TRUE(v1.y() == m.point(vmap[1]).Y);

        EXPECT_TRUE(v2.x() == m.point(vmap[2]).X);
        EXPECT_TRUE(v2.y() == m.point(vmap[2]).Y);

        EXPECT_TRUE(v3.x() == m.point(vmap[3]).X);
        EXPECT_TRUE(v3.y() == m.point(vmap[3]).Y);
    }

    // Test validity, optimality
    {
        edges_are_valid(m);
        edges_are_optimal(m);
    }

    // Test edge and node connections
    {
        EXPECT_TRUE(m.point((size_t)0) == m.point(m.edge(0)));
        EXPECT_TRUE(m.point((size_t)1) == m.point(m.edge(1)));
        EXPECT_TRUE(m.point((size_t)2) == m.point(m.edge(2)));
        EXPECT_TRUE(m.point((size_t)3) == m.point(m.edge(3)));
        EXPECT_TRUE(m.point((size_t)1) == m.point(m.edge(4)));
        EXPECT_TRUE(m.point((size_t)3) == m.point(m.edge(5)));

        for (size_t i = 0; i < 5; i++) {
            const Edge *e = m.edge(i);
            EXPECT_TRUE(e->self() == m.edge(e->next())->prev());
            EXPECT_TRUE(e->self() == m.edge(e->prev())->next());
        }
    }

    // Test triangles
    {
        for (size_t i = 0; i < m.size_triangles(); ++i) {
            Point cc = m.circumcenter(m.triangle(0)->self());
            EXPECT_NEAR(0.5, cc.X, TOL);
            EXPECT_NEAR(0.5, cc.Y, TOL);
        }
    }

    // Forced Refinement
    {
        forced_refinement(m, "square_domain_refine_loop", 7);
    }

    m.delete_me();
}

TEST(Mesh, create__narrow_diamond_domain) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(0.0, 2.0);
    Vertex &v2 = s.new_element<Vertex>(-1.0, 0.0);
    Vertex &v3 = s.new_element<Vertex>(0.0, -2.0);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
    LineSegment &l3 = s.new_element<LineSegment>(v3, v0);

    s.solve();
    s.build();

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, "narrow_diamond_domain");

    // Test number of verticies, edges, triangles
    {
        EXPECT_TRUE(m.size_points() == 4);
        EXPECT_TRUE(m.size_edges() == 6);
        EXPECT_TRUE(m.size_triangles() == 2);

        EXPECT_TRUE(m.num_points() == 4);
        EXPECT_TRUE(m.num_edges() == 5);
        EXPECT_TRUE(m.num_triangles() == 2);
    }

    // Test validity, optimality
    {
        edges_are_valid(m);
        edges_are_optimal(m);
    }

    // Test for proper edge swaps
    {
        std::vector<size_t> vmap = map_verticies_to_points({v0, v1, v2, v3}, m);

        for (size_t i = 5; i < 6; i++) {
            const Edge *e = m.edge(5);

            if (m.point(e->node()) == m.point(vmap[0])) {
                EXPECT_TRUE(m.point(vmap[2]) == m.point(m.edge(e->next())->node()));
                EXPECT_TRUE(m.point(vmap[2]) == m.point(m.edge(e->twin())->node()));
            } else if (m.point(e->node()) == m.point(vmap[2])) {
                EXPECT_TRUE(m.point(vmap[0]) == m.point(m.edge(e->next())->node()));
                EXPECT_TRUE(m.point(vmap[0]) == m.point(m.edge(e->twin())->node()));
            }
        }
    }

    // Test triangle circumcenters
    {
        Point cc0 = m.circumcenter(m.triangle(0)->self()); // TODO: Write m.circumcenter(size_t)
        Point cc1 = m.circumcenter(m.triangle(1)->self());

        EXPECT_NEAR(0.0, cc0.X, TOL);
        EXPECT_NEAR(0.75, std::abs(cc0.Y), TOL);

        EXPECT_NEAR(0.0, cc1.X, TOL);
        EXPECT_NEAR(0.75, std::abs(cc1.Y), TOL);

        EXPECT_NEAR(-cc0.Y, cc1.Y, TOL);
    }

    // Forced Refinement
    {
        forced_refinement(m, "narrow_diamond_domain_refine_loop", 7);
    }

    m.delete_me();
}

TEST(Mesh, create__narrow_rectangle_domain) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(5.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(5.0, 1.0);
    Vertex &v2 = s.new_element<Vertex>(0.0, 2.0);
    Vertex &v3 = s.new_element<Vertex>(-5.0, 1.0);
    Vertex &v4 = s.new_element<Vertex>(-5.0, 0.0);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
    LineSegment &l3 = s.new_element<LineSegment>(v3, v4);
    LineSegment &l4 = s.new_element<LineSegment>(v4, v0);

    s.solve();
    s.build();

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, "narrow_rectangle_domain");

    // Test number of verticies, edges, triangles
    {
        EXPECT_TRUE(m.size_points() == 6);
        EXPECT_TRUE(m.size_edges() == 12);
        EXPECT_TRUE(m.size_triangles() == 4);

        EXPECT_TRUE(m.num_points() == 6);
        EXPECT_TRUE(m.num_edges() == 9);
        EXPECT_TRUE(m.num_triangles() == 4);
    }

    // Test validity, optimality
    {
        edges_are_valid(m);
        edges_are_optimal(m);
    }

    // Test edge splits
    {
        std::vector<size_t> vmap = map_verticies_to_points({v0, v1, v2, v3, v4}, m);

        Point const &v5 = m.point(5);
        EXPECT_NEAR(0.0, v5.X, DBL_EPSILON);
        EXPECT_NEAR(0.0, v5.Y, DBL_EPSILON);

        for (size_t i = 0; i < 12; ++i) {
            Edge const *e = m.edge(i);

            if (m.point(e->node()) == m.point(vmap[2])) {
                if (e->twin() != e->self()) {
                    EXPECT_TRUE(m.point(m.edge(e->twin())->node()) == v5);
                }
            }
        }
    }

    // Forced Refinement
    {
        forced_refinement(m, "narrow_rectangle_domain_refine_loop", 7);
    }

    m.delete_me();
}

TEST(Mesh, create__half_circle_domain) {
    Sketch s;

    Vertex &vc = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v0 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(-1.0, 0.0);

    CircularArc &arc = s.new_element<CircularArc>(v0, v1, vc, 1.0);
    LineSegment &line = s.new_element<LineSegment>(v1, v0);

    s.solve();
    s.build();

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, "half_circle_domain");

    // Test number of verticies, edges, triangles
    {
        EXPECT_TRUE(m.size_points() == 6);
        EXPECT_TRUE(m.size_edges() == 12);
        EXPECT_TRUE(m.size_triangles() == 4);

        EXPECT_TRUE(m.num_points() == 6);
        EXPECT_TRUE(m.num_edges() == 9);
        EXPECT_TRUE(m.num_triangles() == 4);
    }

    // Test validity, optimality
    {
        edges_are_valid(m);
        edges_are_optimal(m);
    }

    // Forced Refinement
    {
        forced_refinement(m, "half_circle_domain_refine_loop", 7);
    }

    m.delete_me();
}

TEST(Mesh, create__horseshoe_domain) {
    Sketch s;

    Vertex &vc = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v0 = s.new_element<Vertex>(2.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(-2.0, 0.0);
    Vertex &v2 = s.new_element<Vertex>(-1.9, 0.0);
    Vertex &v3 = s.new_element<Vertex>(1.9, 0.0);

    CircularArc &arc0 = s.new_element<CircularArc>(v0, v1, vc, 2.0);
    LineSegment &line0 = s.new_element<LineSegment>(v1, v2);
    CircularArc &arc1 = s.new_element<CircularArc>(v3, v2, vc, 1.9);
    LineSegment &line1 = s.new_element<LineSegment>(v3, v0);

    s.solve();
    s.build();

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, "horseshoe_domain");

    { // Test triangles, possibly redundant
        for (size_t i = 0; i < m.size_edges(); ++i) {
            EXPECT_TRUE(m.edge(m.edge(i)->next())->prev() == m.edge(i)->self());
            EXPECT_TRUE(m.edge(m.edge(m.edge(m.edge(i)->next())->next())->prev())->prev() == m.edge(i)->self());
            EXPECT_TRUE(m.edge(m.edge(m.edge(i)->next())->next())->next() == m.edge(i)->self());

            EXPECT_TRUE(m.edge(m.edge(i)->prev())->next() == m.edge(i)->self());
            EXPECT_TRUE(m.edge(m.edge(m.edge(m.edge(i)->prev())->prev())->next())->next() == m.edge(i)->self());
            EXPECT_TRUE(m.edge(m.edge(m.edge(i)->prev())->prev())->prev() == m.edge(i)->self());

            if ((m.edge(i)->twin() != m.edge(i)->self())) {
                EXPECT_TRUE(m.edge(m.edge(i)->twin())->node() == m.edge(m.edge(i)->next())->node());
            }
        }
    }

    // Test validity, optimality
    {
        edges_are_valid(m);
        edges_are_optimal(m);
    }

    {
        forced_refinement(m, "horseshoe_domain_refine_loop", 7);
    }

    m.delete_me();
}

TEST(Mesh, create__I_shaped_domain) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(3.0, 0.0);
    Vertex &v2 = s.new_element<Vertex>(3.0, 1.0);
    Vertex &v3 = s.new_element<Vertex>(2.0, 1.0);
    Vertex &v4 = s.new_element<Vertex>(2.0, 2.0);
    Vertex &v5 = s.new_element<Vertex>(3.0, 2.0);
    Vertex &v6 = s.new_element<Vertex>(3.0, 3.0);
    Vertex &v7 = s.new_element<Vertex>(0.0, 3.0);
    Vertex &v8 = s.new_element<Vertex>(0.0, 2.0);
    Vertex &v9 = s.new_element<Vertex>(1.0, 2.0);
    Vertex &v10 = s.new_element<Vertex>(1.0, 1.0);
    Vertex &v11 = s.new_element<Vertex>(0.0, 1.0);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
    LineSegment &l3 = s.new_element<LineSegment>(v3, v4);
    LineSegment &l4 = s.new_element<LineSegment>(v4, v5);
    LineSegment &l5 = s.new_element<LineSegment>(v5, v6);
    LineSegment &l6 = s.new_element<LineSegment>(v6, v7);
    LineSegment &l7 = s.new_element<LineSegment>(v7, v8);
    LineSegment &l8 = s.new_element<LineSegment>(v8, v9);
    LineSegment &l9 = s.new_element<LineSegment>(v9, v10);
    LineSegment &l10 = s.new_element<LineSegment>(v10, v11);
    LineSegment &l11 = s.new_element<LineSegment>(v11, v0);

    s.solve();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "i_shaped_domain");

    s.build();

    EXPECT_TRUE(s.size_contours() == 1);

    Mesh m{s};

    m.create();

    m.save_as(SAVE_DIR, "i_shaped_domain_mesh");

    // Test validity, optimality
    {
        edges_are_valid(m);
        edges_are_optimal(m);
    }

    // Forced Refinement
    {
        forced_refinement(m, "i_shaped_domain_mesh_refine_loop", 7);
    }

    m.delete_me();
}

TEST(Mesh, create__corner_square_domain) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(0.5, 0.0);
    Vertex &v2 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v3 = s.new_element<Vertex>(1.0, 0.5);
    Vertex &v4 = s.new_element<Vertex>(1.0, 1.0);
    Vertex &v5 = s.new_element<Vertex>(0.0, 1.0);
    Vertex &v6 = s.new_element<Vertex>(0.5, 0.5);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v6);
    LineSegment &l2 = s.new_element<LineSegment>(v6, v3);
    LineSegment &l3 = s.new_element<LineSegment>(v3, v4);
    LineSegment &l4 = s.new_element<LineSegment>(v4, v5);
    LineSegment &l5 = s.new_element<LineSegment>(v5, v0);

    LineSegment &l6 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l7 = s.new_element<LineSegment>(v2, v3);

    s.solve();
    s.build();

    Mesh m{s};
    m.create();

    m.save_as(SAVE_DIR, "corner_square_domain");

    // Test number of vertices, edges, triangles
    {
        EXPECT_TRUE(m.num_points() == 9);
        EXPECT_TRUE(m.size_points() == 9);

        EXPECT_TRUE(m.num_edges() == 16);
        EXPECT_TRUE(m.size_edges() == 24);

        EXPECT_TRUE(m.num_triangles() == 8);
        EXPECT_TRUE(m.size_triangles() == 8);
    }

    { // Test validity, optimality
        edges_are_valid(m);
        edges_are_optimal(m);
    }

    { // Forced refinement
        forced_refinement(m, "corner_square_domain_refine_loop", 7);
    }

    m.delete_me();
}

TEST(Mesh, create__square_in_square_domain) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(2.0, 0.0);
    Vertex &v2 = s.new_element<Vertex>(2.0, 2.0);
    Vertex &v3 = s.new_element<Vertex>(0.0, 2.0);

    Vertex &v4 = s.new_element<Vertex>(0.9, 0.9);
    Vertex &v5 = s.new_element<Vertex>(1.9, 0.9);
    Vertex &v6 = s.new_element<Vertex>(1.9, 1.9);
    Vertex &v7 = s.new_element<Vertex>(0.9, 1.9);

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
    LineSegment &l3 = s.new_element<LineSegment>(v3, v0);

    LineSegment &l4 = s.new_element<LineSegment>(v4, v5);
    LineSegment &l5 = s.new_element<LineSegment>(v5, v6);
    LineSegment &l6 = s.new_element<LineSegment>(v6, v7);
    LineSegment &l7 = s.new_element<LineSegment>(v7, v4);

    s.solve();
    s.build();

    s.save_as<SaveMethod::Rasterize>(SAVE_DIR, "square_in_square_domain");

    Mesh m{s};
    m.MaximumElementSize = 0.1;
    m.MinimumElementSize = 0.01;
    m.MinimumElementQuality = M_SQRT1_2;

    m.create();

    m.save_as(SAVE_DIR, "square_in_square_domain_mesh");

    edges_are_valid(m);
    edges_are_optimal(m);

    m.refine();

    m.save_as(SAVE_DIR, "square_in_square_domain_mesh_refine");

    m.delete_me();
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

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(2.0, 0.0);
    Vertex &v2 = s.new_element<Vertex>(1.0, sqrt(3.0));

    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l2 = s.new_element<LineSegment>(v2, v0);

    s.solve();
    EXPECT_TRUE(s.build());

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
        Edge const *e;
        Point vp;
        size_t loc = i;

        vp = vi;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Interior);

        vp = ve0;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[0])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[1])));

        vp = ve1;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[1])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[2])));

        vp = ve2;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[2])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[0])));
    }

    // Corner Cases
    ve0 = {1.0, -DBL_EPSILON};
    ve1 = {1.5 + DBL_EPSILON, sqrt(3.0) / 2.0 + DBL_EPSILON};
    ve2 = {0.5 - DBL_EPSILON, sqrt(3.0) / 2.0 + DBL_EPSILON};

    for (size_t i = 0; i != mesh.size_edges(); ++i) {
        Edge const *e;
        Point vp;
        size_t loc = i;

        vp = ve0;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[0])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[1])));

        vp = ve1;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[1])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[2])));

        vp = ve2;
        EXPECT_TRUE(mesh.locate_triangle(vp, loc) == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[2])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[0])));
    }
}

TEST(Mesh, locate_triange__square_domain) {
    Sketch s;

    Vertex &v0 = s.new_element<Vertex>(0.0, 0.0);
    Vertex &v1 = s.new_element<Vertex>(1.0, 0.0);
    Vertex &v2 = s.new_element<Vertex>(1.0, 1.0);
    Vertex &v3 = s.new_element<Vertex>(0.0, 1.0);
    LineSegment &l0 = s.new_element<LineSegment>(v0, v1);
    LineSegment &l1 = s.new_element<LineSegment>(v1, v2);
    LineSegment &l2 = s.new_element<LineSegment>(v2, v3);
    LineSegment &l3 = s.new_element<LineSegment>(v3, v0);

    s.solve();
    s.build();
    /*
    std::vector<const Curve*> cc{ &l0,&l1,&l2,&l3 };
    Contour cont{ cc };

    std::vector<const Contour*> cv{ &cont };
    Mesh mesh{ cv };
    */
    Mesh mesh{s};
    mesh.create();

    EXPECT_TRUE(mesh.size_edges() == 6);
    EXPECT_TRUE(mesh.size_points() == 4);

    // Interior Points
    Point vi0{0.25, 0.25};
    Point vi1{0.75, 0.75};
    for (size_t i = 0; i < mesh.size_edges(); ++i) {
        //Edge const *e;
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
        Edge const *e;
        Point vp;
        size_t loc = i;

        vp = ve0;
        EXPECT_TRUE(mesh.locate_triangle(vp,loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[0])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[1])));

        vp = ve1;
        EXPECT_TRUE(mesh.locate_triangle(vp,loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[1])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[2])));

        vp = ve2;
        EXPECT_TRUE(mesh.locate_triangle(vp,loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[2])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[3])));

        vp = ve3;
        EXPECT_TRUE(mesh.locate_triangle(vp,loc) == LocateTriangleResult::Exterior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[3])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[0])));
    }

    // Exterior Edge Points
    ve0 = {0.5, -DBL_EPSILON};
    ve1 = {1.0 + DBL_EPSILON, 0.5};
    ve2 = {0.5, 1.0 + DBL_EPSILON};
    ve3 = {0.0 - DBL_EPSILON, 0.5};
    LocateTriangleResult result;
    for (size_t i = 0; i < mesh.size_edges(); ++i) {
        Edge const *e;
        Point vp;
        size_t loc = i;

        vp = ve0;
        result = mesh.locate_triangle(vp,loc);
        EXPECT_TRUE(result == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[0])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[1])));

        vp = ve1;
        result = mesh.locate_triangle(vp,loc);
        EXPECT_TRUE(result == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[1])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[2])));

        vp = ve2;
        result = mesh.locate_triangle(vp,loc);
        EXPECT_TRUE(result == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[2])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[3])));

        vp = ve3;
        result = mesh.locate_triangle(vp,loc);
        EXPECT_TRUE(result == LocateTriangleResult::Interior);
        e = mesh.edge(loc);
        EXPECT_TRUE((mesh.point(e->node()) == mesh.point(vmap[3])) && (mesh.point(mesh.edge(e->next())->node()) == mesh.point(vmap[0])));
    }

    Point vie{0.5, 0.5};
    for (size_t i = 0; i < mesh.size_edges(); ++i) {
        //const Edge *e = mesh.edge(i);

        EXPECT_TRUE(mesh.locate_triangle(vie, i) == LocateTriangleResult::Interior);
    }
}