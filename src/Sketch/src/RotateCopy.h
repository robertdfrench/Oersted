#ifndef OERSTED_ROTATECOPY_H
#define OERSTED_ROTATECOPY_H

#include "Pattern.h"

class RotateCopy : public Pattern {
public:
    RotateCopy(std::vector<const Curve *> &input, Vertex *center, double angle, size_t copies);

private:
    Vertex *Center;
    double Angle;
    size_t Copies;
};

#endif //OERSTED_ROTATECOPY_H