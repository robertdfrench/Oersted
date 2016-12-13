#include "Rotation.h"
#include "Vertex.h"

void Rotation::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
    double xo = Origin->x();
    double yo = Origin->y();

    double x0 = V0->x();
    double y0 = V0->y();

    double x1 = V1->x();
    double y1 = V1->y();

    double a = Angle * M_PI / 180.0;

    size_t ei = EquationIndex;
    r(ei) = (x0 - xo) * cos(a) - (y0 - yo) * sin(a) - (x1 - xo);

    J(ei, Origin->x_index()) = 1.0 - cos(a);
    J(ei, Origin->y_index()) = sin(a);
    J(ei, V0->x_index()) = cos(a);
    J(ei, V0->y_index()) = -sin(a);
    J(ei, V1->x_index()) = -1.0;

    ei++;
    r(ei) = (x0 - xo) * sin(a) + (y0 - yo) * cos(a) - (y1 - yo);

    J(ei, Origin->x_index()) = -sin(a);
    J(ei, Origin->y_index()) = 1 - cos(a);
    J(ei, V0->x_index()) = sin(a);
    J(ei, V0->y_index()) = cos(a);
    J(ei, V1->y_index()) = -1.0;
}