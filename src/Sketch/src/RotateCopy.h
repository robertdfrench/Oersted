#ifndef OERSTED_ROTATECOPY_H
#define OERSTED_ROTATECOPY_H

#include "Pattern.h"

class RotateCopy : public Pattern {
public:
    RotateCopy(std::vector<std::shared_ptr<Curve>> &input, std::shared_ptr<Vertex> center, double angle, size_t copies, bool remove_internal = false);

private:
    std::shared_ptr<Vertex> Center;
    double Angle; // TODO: Write 'update_angle' method which updates all associated rotation constraints
    size_t Copies;
};

#endif //OERSTED_ROTATECOPY_H