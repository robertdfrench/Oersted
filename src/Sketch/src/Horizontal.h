#ifndef OERSTED_HORIZONTAL_H
#define OERSTED_HORIZONTAL_H

#include "Constraint.h"

class LineSegment;

class Horizontal : public Constraint {
public:
    Horizontal(std::shared_ptr<LineSegment const> l) : Line(l) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override;

protected:
    std::shared_ptr<LineSegment const> Line;
};

#endif //OERSTED_HORIZONTAL_H
