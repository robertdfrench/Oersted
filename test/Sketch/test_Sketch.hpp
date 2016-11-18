#ifndef OERSTED_TEST_SKETCH_H
#define OERSTED_TEST_SKETCH_H

#include "Sketch.hpp"
#include "gtest.h"

#include <cmath>

#define TOL FLT_EPSILON //#TODO: Tolerance is limited by accuracy of tangency constraint
#define SAVE_DIR "./test/output/Sketch/"

void test_sketch_size(Sketch &s, size_t nverts, size_t ncurves, size_t nconstraints, size_t ncontours);

#endif