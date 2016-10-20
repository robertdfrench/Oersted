#include "Sketch.hpp"

Fixation::Fixation(Vertex &v) {
    Point = &v;
    Dim = new Vertex(v);
}

void Fixation::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    r(EquationIndex) = Point->x() - Dim->x();
    J(EquationIndex, Point->X->get_index()) += 1.0;

    r(EquationIndex + 1) = Point->y() - Dim->y();
    J(EquationIndex + 1, Point->Y->get_index()) += 1.0;
}