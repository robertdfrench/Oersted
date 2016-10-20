#include "Sketch.hpp"

void Horizontal::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    r(EquationIndex) = Line->Start->y() - Line->End->y();

    J(EquationIndex, Line->Start->Y->get_index()) += 1.0;
    J(EquationIndex, Line->End->Y->get_index()) -= 1.0;
}