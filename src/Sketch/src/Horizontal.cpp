#include "Horizontal.h"
#include "LineSegment.h"
#include "Vertex.h"

void Horizontal::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
    r(EquationIndex) = Line->start()->y() - Line->end()->y();

    J(EquationIndex, Line->start()->y_index()) += 1.0;
    J(EquationIndex, Line->end()->y_index()) -= 1.0;
}