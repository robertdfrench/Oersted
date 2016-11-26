#ifndef OERSTED_PATTERN_H
#define OERSTED_PATTERN_H

#include "SketchElement.h"

class Curve;
class Constraint;
class Vertex;

class Pattern : public SketchElement {
public:
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 0;
    };

    void register_elements(Sketch *s);

    void register_parameters(Sketch *s) override {};

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override {};

protected:
    std::vector<std::shared_ptr<Curve>> Input;

    std::vector<std::shared_ptr<Constraint>> Constraints;
    std::vector<std::shared_ptr<Curve>> Curves;
    std::vector<std::shared_ptr<Vertex>> Verticies;

    bool RemoveInternalBoundaries;
};

#endif //OERSTED_PATTERN_H