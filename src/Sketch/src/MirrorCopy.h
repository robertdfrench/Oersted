#ifndef OERSTED_MIRRORCOPY_H
#define OERSTED_MIRRORCOPY_H

#include "Pattern.h"

class LineSegment;

class MirrorCopy : public Pattern {
public:
    MirrorCopy(std::vector<std::shared_ptr<Curve const>> &input, std::shared_ptr<LineSegment const> l, bool remove_internal = false);

protected:
    std::shared_ptr<LineSegment const> SymmetryLine;
};

#endif //OERSTED_MIRRORCOPY_H