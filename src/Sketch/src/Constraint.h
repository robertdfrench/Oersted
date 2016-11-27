#ifndef OERSTED_CONSTRAINT_H
#define OERSTED_CONSTRAINT_H

#include "SketchElement.h"

class Constraint : public SketchElement {
public:
    void register_parameters(Sketch *s) const override {};
};

#endif //OERSTED_CONSTRAINT_H