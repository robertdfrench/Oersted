#include "test_Sketch.hpp"

TEST(Vertex, constructor) {
    { // ARGS::()
        EXPECT_NO_THROW(Vertex v);
    }

    { // ARGS::(double,double)
        EXPECT_NO_THROW(Vertex v(3.14159, 2.71828));
    }
}

TEST(Vertex, rotate) {
    std::shared_ptr<Vertex> v = std::make_shared<Vertex>(1.0, 1.0);
    std::shared_ptr<Vertex> origin = std::make_shared<Vertex>(2.0, 2.0);
    double angle{45.0};
    double r{hypot(v->x() - origin->x(), v->y() - origin->y())};

    double x, y;

    std::tie(x, y) = v->rotate(origin, angle);

    EXPECT_NEAR(hypot(x - origin->x(), y - origin->y()), r, r * TOL);
    EXPECT_NEAR(origin->x(), x, r * TOL);
    EXPECT_NEAR(origin->y() - r, y, r * TOL);
}