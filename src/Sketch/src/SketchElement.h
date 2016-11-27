#ifndef OERSTED_SKETCHELEMENT_H
#define OERSTED_SKETCHELEMENT_H

#define SIGN(x) (double)((x > 0.0) - (x < 0.0))

#include <cmath>
#include <memory>
#include <utility>
#include <vector>
#include <list>

#include "Eigen"

class Sketch;

class SketchElement {
public:
    size_t get_equation_index() const { return EquationIndex; };

    virtual size_t set_equation_index(size_t i) = 0;

    virtual void register_parameters(Sketch *s) const = 0;

    virtual void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const = 0;

protected:
    size_t EquationIndex;
};

#endif //OERSTED_SKETCHELEMENT_H
