#ifndef OERSTED_SYMMETRY_H
#define OERSTED_SYMMETRY_H

#include "Constraint.h"

class Vertex;
class LineSegment;

class Symmetry : public Constraint {
public:
    Symmetry(std::shared_ptr<Vertex const> v0, std::shared_ptr<Vertex const> v1, std::shared_ptr<LineSegment const> line) : V0(v0), V1(v1), SymmetryLine(line) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override;

protected:
    std::shared_ptr<Vertex const> V0;

    std::shared_ptr<Vertex const> V1;

    std::shared_ptr<LineSegment const> SymmetryLine;
};

#endif //OERSTED_SYMMETRY_H
