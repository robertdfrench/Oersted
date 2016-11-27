#include "Radius.h"
#include "CircularArc.h"

void Radius::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
    r(EquationIndex) = Arc->radius() - Dim;

    J(EquationIndex, Arc->radius_index()) += 1.0;
}