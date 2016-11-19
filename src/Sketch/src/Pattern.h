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
    std::vector<std::shared_ptr<Curve>> Input;
    bool RemoveInternalBoundaries;

    std::vector<std::shared_ptr<Vertex>> Verticies;
    std::vector<std::shared_ptr<Curve>> Curves;
    std::vector<std::shared_ptr<Constraint>> Constraints;
};

#include "MirrorCopy.h"
#include "RotateCopy.h"

#endif //OERSTED_PATTERN_H