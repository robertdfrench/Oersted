#ifndef OERSTED_ANGLE_H
#define OERSTED_ANGLE_H

#include "Constraint.h"

class LineSegment;

class Angle : public Constraint {
public:
    Angle(std::shared_ptr<LineSegment> l0, std::shared_ptr<LineSegment> l1, double angle) : Line0(l0), Line1(l1), Dim(angle) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    double dim() const { return Dim; };

    void dim(double d) { Dim = d;};

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;

protected:
    std::shared_ptr<LineSegment> Line0;
    std::shared_ptr<LineSegment> Line1;

    double Dim;
};

#endif //OERSTED_ANGLE_H