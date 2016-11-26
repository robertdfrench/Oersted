#ifndef OERSTED_FIXATION_H
#define OERSTED_FIXATION_H

#include "Constraint.h"
#include "sPoint.h"

class Vertex;

class Fixation : public Constraint {
public:
    Fixation(std::shared_ptr<Vertex> v);

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;

protected:
    std::shared_ptr<Vertex> Point;
    sPoint Dim;
};

#endif //OERSTED_FIXATION_H