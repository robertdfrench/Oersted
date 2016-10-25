#ifndef OERSTED_UTIL_H
#define OERSTED_UTIL_H

#include "Mesh.h"

bool are_intersecting(const Edge *e0, const Edge *e1);

void circumradius(std::vector<Edge *> &triangles, std::vector<double> &radii);

void shortest_edge_length(std::vector<Edge *> &triangles, std::vector<double> &length);

void element_quality(std::vector<Edge *> &triangles, std::vector<double> &radii, std::vector<double> &quality);

void sort_permutation(std::vector<double> &value, std::vector<size_t> &index);

#endif //OERSTED_UTIL_H
