#ifndef OERSTED_FIXATION_H
#define OERSTED_FIXATION_H

#include "Constraint.h"
#include "doublen.h"

class Vertex;

class Fixation : public Constraint {
public:
    Fixation(std::shared_ptr<Vertex const> v);

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override;

protected:
    double2 Dim;
    std::shared_ptr<Vertex const> Point;
};

#endif //OERSTED_FIXATION_H