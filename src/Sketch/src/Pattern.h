#ifndef OERSTED_PATTERN_H
#define OERSTED_PATTERN_H

#include "Sketch.h"

class Pattern : public SketchElement {
public:
    void register_elements(Sketch *s);

    void register_parameters(Sketch *s) override {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 0;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override {};

protected:
    std::vector<const Curve *> Input;
    bool RemoveInternalBoundaries;

    std::vector<Vertex *> Verticies;
    std::vector<Curve *> Curves;
    std::vector<Constraint *> Constraints;
};

#include "MirrorCopy.h"
#include "RotateCopy.h"

#endif //OERSTED_PATTERN_H