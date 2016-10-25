#ifndef OERSTED_POINT_H
#define OERSTED_POINT_H

#include "Mesh.h"

class Point { // #TODO, Replace Verticies member in Mesh with lighterweight Node class
public:

    Point() : X{0.0}, Y{0.0} {};

    Point(double x, double y) : X{x}, Y{y} {};

    Point(const Vertex &v) : X{v.x()}, Y{v.y()} {};

    Point(const Vertex *v) : X{v->x()}, Y{v->y()} {};

    //double W; // Nurbs weight?
    double X;
    double Y;
    //double Z; // 3-Dimensions?

    //Vertex* V = nullptr; Keep track of constrained nodes?

    bool operator==(const Point &p) const { return (X == p.X) && (Y == p.Y); };

    bool operator==(const Vertex &v) const { return (X == v.x()) && (Y == v.y()); };

    bool operator!=(const Point &p) const { return (X != p.X) || (Y != p.Y); };

    bool operator!=(const Vertex &v) const { return (X != v.x()) && (Y != v.y()); };
};

#endif //OERSTED_POINT_H
