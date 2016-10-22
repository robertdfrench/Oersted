#ifndef OERSTED_CONSTRAINT_H
#define OERSTED_CONSTRAINT_H

#include "Sketch.h"

class Constraint : public SketchElement {
public:
    void register_parameters(Sketch *s) override {};
};

#include "Angle.h"
#include "Coincident.h"
#include "Distance.h"
#include "Constraint.h"
#include "Fixation.h"
#include "Horizontal.h"
#include "Length.h"
#include "Radius.h"
#include "Rotation.h"
#include "Symmetry.h"
#include "Tangency.h"
#include "Vertical.h"

#endif //OERSTED_CONSTRAINT_H