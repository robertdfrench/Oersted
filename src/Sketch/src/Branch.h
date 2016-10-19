#ifndef OERSTED_BRANCH_H
#define OERSTED_BRANCH_H

struct Branch {
    const Curve *Path;
    bool Orientation;
    double Angle;
};

#endif //OERSTED_BRANCH_H
