#ifndef OERSTED_ROTATECOPY_H
#define OERSTED_ROTATECOPY_H

#include "Pattern.h"

class RotateCopy : public Pattern {
public:
    RotateCopy(std::vector<std::shared_ptr<Curve const>> input, std::shared_ptr<Vertex const> center, double angle, size_t copies, bool remove_internal = false);

    double angle() const { return Angle; };

    //void angle(double a) { Angle = a; }; // TODO: Must update all associated rotation constraints

private:
    std::shared_ptr<Vertex const> Center;

    double Angle;

    size_t Copies;
};

#endif //OERSTED_ROTATECOPY_H