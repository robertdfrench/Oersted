#ifndef OERSTED_VERTICAL_H
#define OERSTED_VERTICAL_H

#include "Constraint.h"

class LineSegment;

class Vertical : public Constraint {
public:
    Vertical(std::shared_ptr<LineSegment const> l) : Line(l) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override;

protected:
    std::shared_ptr<LineSegment const> Line;
};

#endif //OERSTED_VERTICAL_H
