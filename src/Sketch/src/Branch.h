#ifndef OERSTED_BRANCH_H
#define OERSTED_BRANCH_H

#include <memory>

class Curve;

struct Branch {
    std::shared_ptr<Curve const> Path;

    double Angle;

    bool Orientation;
};

#endif //OERSTED_BRANCH_H
