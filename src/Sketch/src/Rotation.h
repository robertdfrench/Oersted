#ifndef OERSTED_ROTATION_H
#define OERSTED_ROTATION_H

#include "Constraint.h"

class Vertex;

class Rotation : public Constraint {
public:
    Rotation(std::shared_ptr<Vertex const> v0, std::shared_ptr<Vertex const> v1, std::shared_ptr<Vertex const> origin, double a) : V0(v0), V1(v1), Origin(origin), Angle(a) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override;

protected:
    std::shared_ptr<Vertex const> V0;

    std::shared_ptr<Vertex const> V1;

    std::shared_ptr<Vertex const> Origin;

    double Angle;
};

#endif //OERSTED_ROTATION_H
