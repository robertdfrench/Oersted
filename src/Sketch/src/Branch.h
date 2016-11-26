#ifndef OERSTED_BRANCH_H
#define OERSTED_BRANCH_H

#include <memory>

class Curve;

struct Branch {
    std::shared_ptr<Curve> Path;
    double Angle;
    bool Orientation;
};

#endif //OERSTED_BRANCH_H
