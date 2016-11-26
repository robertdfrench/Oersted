#ifndef OERSTED_SPOINT_H
#define OERSTED_SPOINT_H

#include <cfloat>

class sPoint {
public:
    sPoint() : X{DBL_MAX}, Y{DBL_MAX} {};
    sPoint(double x, double y) : X{x}, Y{y} {};

    double X;
    double Y;

    double x() const {return X;};
    double y() const {return Y;};
};

#endif
