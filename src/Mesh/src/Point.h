#ifndef OERSTED_POINT_H
#define OERSTED_POINT_H

#include "Mesh.h"

class Point { // TODO, Replace Verticies member in Mesh with lighterweight Node class
public:

    Point() : X{0.0}, Y{0.0} {};

    Point(double x, double y) : X{x}, Y{y} {};

    Point(Vertex const &v) : X{v.x()}, Y{v.y()} {};

    Point(Vertex const *v) : X{v->x()}, Y{v->y()} {};

    //double W; // Nurbs weight?
    double X;
    double Y;
    //double Z; // 3-Dimensions?

    bool operator==(Point const &p) const { return (X == p.X) && (Y == p.Y); };

    bool operator==(Vertex const &v) const { return (X == v.x()) && (Y == v.y()); };

    bool operator!=(Point const &p) const { return (X != p.X) || (Y != p.Y); };

    bool operator!=(Vertex const &v) const { return (X != v.x()) && (Y != v.y()); };
};

double dist(Point const &p0, Point const &p1);

#endif //OERSTED_POINT_H
