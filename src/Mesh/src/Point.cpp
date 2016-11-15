#include "Mesh.hpp"

double dist(Point const &p0, Point const &p1) {
    return std::sqrt((p0.X - p1.X) * (p0.X - p1.X) + (p0.Y - p1.Y) * (p0.Y - p1.Y));
};