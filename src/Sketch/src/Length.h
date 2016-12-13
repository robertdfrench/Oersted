#ifndef OERSTED_LENGTH_H
#define OERSTED_LENGTH_H

#include "Constraint.h"

class LineSegment;

class Length : public Constraint {
public:
    Length(std::shared_ptr<LineSegment const> c, double length) : Line(c), Dim(length) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override;

protected:
    std::shared_ptr<LineSegment const> Line;

    double Dim;
};

#endif //OERSTED_LENGTH_H