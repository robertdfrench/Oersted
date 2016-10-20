#include "Sketch.hpp"

void Vertical::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    r(EquationIndex) = Line->Start->x() - Line->End->x();

    J(EquationIndex, Line->Start->X->get_index()) += 1.0;
    J(EquationIndex, Line->End->X->get_index()) -= 1.0;
}