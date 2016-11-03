#ifndef OERSTED_MESH_TEST_UTIL_H
#define OERSTED_MESH_TEST_UTIL_H

bool edges_are_optimal(Mesh &m);
bool edges_are_valid(Mesh &m);
void forced_refinement(Mesh &m, std::string file_name, size_t num_refines);
std::vector<size_t> map_verticies_to_points(std::vector<Vertex> verts, Mesh m);

#endif //OERSTED_UTIL_H
