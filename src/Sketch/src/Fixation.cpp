#include "Sketch.hpp"

Fixation::Fixation(std::shared_ptr<Vertex> v) {
    Point = v;
    Dim = sPoint(v->x(), v->y());
}

void Fixation::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    r(EquationIndex) = Point->x() - Dim.X;
    J(EquationIndex, Point->X->get_index()) += 1.0;

    r(EquationIndex + 1) = Point->y() - Dim.Y;
    J(EquationIndex + 1, Point->Y->get_index()) += 1.0;
}