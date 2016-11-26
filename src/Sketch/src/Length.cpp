#include "Length.h"
#include "LineSegment.h"
#include "Vertex.h"

void Length::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    double dx = Line->end()->x() - Line->start()->x();
    double dy = Line->end()->y() - Line->start()->y();
    const double d = sqrt(dx * dx + dy * dy);

    r(EquationIndex) = d - Dim;

    dx = dx / d;
    dy = dy / d;

    J(EquationIndex, Line->start()->x_index()) -= dx;
    J(EquationIndex, Line->start()->y_index()) -= dy;

    J(EquationIndex, Line->end()->x_index()) += dx;
    J(EquationIndex, Line->end()->y_index()) += dy;
}