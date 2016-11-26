#include "Vertical.h"
#include "LineSegment.h"

void Vertical::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    r(EquationIndex) = Line->start()->x() - Line->end()->x();

    J(EquationIndex, Line->start()->x_index()) += 1.0;
    J(EquationIndex, Line->end()->x_index()) -= 1.0;
}