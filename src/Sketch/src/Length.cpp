#include "Sketch.hpp"

void Length::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    double dx = Line->End->x() - Line->Start->x();
    double dy = Line->End->y() - Line->Start->y();
    const double d = sqrt(dx * dx + dy * dy);

    r(EquationIndex) = d - Dim;

    dx = dx / d;
    dy = dy / d;

    J(EquationIndex, Line->Start->X->get_index()) -= dx;
    J(EquationIndex, Line->Start->Y->get_index()) -= dy;

    J(EquationIndex, Line->End->X->get_index()) += dx;
    J(EquationIndex, Line->End->Y->get_index()) += dy;
}