#ifndef OERSTED_RADIUS_H
#define OERSTED_RADIUS_H

#include "Constraint.h"

class CircularArc;

class Radius : public Constraint {
public:
    Radius(std::shared_ptr<CircularArc> c, double r) : Arc(c), Dim(r) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    double dim() const { return Dim; };

    void dim(double d) { Dim = d; };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;

protected:
    std::shared_ptr<CircularArc> Arc;

    double Dim;
};

#endif //OERSTED_RADIUS_H
