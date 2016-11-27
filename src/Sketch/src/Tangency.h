#ifndef OERSTED_TANGENCY_H
#define OERSTED_TANGENCY_H

#include "Constraint.h"

class CircularArc;
class LineSegment;

class Tangency : public Constraint { // TODO: Template tangency for CircularArc/CircularArc and CircularArc/LineSegment tangency
public:
    Tangency(std::shared_ptr<CircularArc const> ca, std::shared_ptr<LineSegment const> ls) : Arc(ca), Line(ls) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override;

protected:
    std::shared_ptr<CircularArc const> Arc;

    std::shared_ptr<LineSegment const> Line;
};

#endif //OERSTED_TANGENCY_H
