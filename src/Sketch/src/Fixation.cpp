#include "Fixation.h"
#include "Vertex.h"

Fixation::Fixation(std::shared_ptr<Vertex const> v) {
    Point = v;
    Dim = double2{v->x(), v->y()};
}

void Fixation::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
    r(EquationIndex) = Point->x() - Dim.X;
    J(EquationIndex, Point->x_index()) += 1.0;

    r(EquationIndex + 1) = Point->y() - Dim.Y;
    J(EquationIndex + 1, Point->y_index()) += 1.0;
}