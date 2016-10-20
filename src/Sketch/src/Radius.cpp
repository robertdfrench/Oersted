#include "Sketch.hpp"

void Radius::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    r(EquationIndex) = Arc->radius() - Dim;

    J(EquationIndex, Arc->Radius->get_index()) += 1.0;
}