#ifndef OERSTED_MESH_UTIL_H
#define OERSTED_MESH_UTIL_H

#include "Mesh.h"

// TODO: Make these mesh routines
bool are_intersecting(Edge const *e0, Edge const *e1, Mesh const &m);

bool in_triangle(Point const *p, Edge const *&e, Mesh const &m);

void element_quality(std::vector<Edge *> &triangles, std::vector<double> &radii, std::vector<double> &quality, Mesh const &m);

void sort_permutation_ascending(std::vector<double> &value, std::vector<size_t> &index);

void sort_permutation_descending(std::vector<double> &values, std::vector<size_t> &index);

#endif //OERSTED_MESH_UTIL_H
