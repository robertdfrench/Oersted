#ifndef OERSTED_MIRRORCOPY_H
#define OERSTED_MIRRORCOPY_H

#include "Pattern.h"

class MirrorCopy : public Pattern {
public:
    MirrorCopy(std::vector<const Curve *> &input, LineSegment *l, bool remove_internal = false);

private:
    LineSegment *SymmetryLine;
};

#endif //OERSTED_MIRRORCOPY_H