#ifndef OERSTED_MIRRORCOPY_H
#define OERSTED_MIRRORCOPY_H

#include "Pattern.h"

class MirrorCopy : public Pattern {
public:
    MirrorCopy(std::vector<std::shared_ptr<Curve>> &input, std::shared_ptr<LineSegment> l, bool remove_internal = false);

private:
    std::shared_ptr<LineSegment> SymmetryLine;
};

#endif //OERSTED_MIRRORCOPY_H