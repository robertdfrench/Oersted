#ifndef OERSTED_TEST_SKETCH_H
#define OERSTED_TEST_SKETCH_H

#include <cmath>
#include <ctgmath>

#include "Sketch.hpp"
#include "gtest.h"

#define TOL FLT_EPSILON // TODO: Tolerance is limited by accuracy of tangency constraint
#define SAVE_DIR "./test/output/Sketch/"

void test_sketch_size(Sketch &s, size_t nverts, size_t ncurves, size_t nconstraints, size_t ncontours);

bool has_rotational_image(Sketch &s, std::shared_ptr<Vertex const> v, std::shared_ptr<Vertex const> center, double angle);

bool has_rotational_image(Sketch &s, std::shared_ptr<Vertex const> v0, std::shared_ptr<Vertex const> v1, std::shared_ptr<Vertex const> center, double angle);

void test_rotated_verticies(Sketch &s, std::vector<size_t> index, std::shared_ptr<Vertex const> center, double angle, size_t copies);

void test_rotated_curves(Sketch &s, std::vector<size_t> index, std::shared_ptr<Vertex const> center, double angle, size_t copies);

#endif