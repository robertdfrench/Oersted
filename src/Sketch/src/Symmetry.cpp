#include "Sketch.hpp"

void Symmetry::update(Eigen::MatrixXd &J, Eigen::VectorXd &r) {
    double x0 = V0->x();
    double y0 = V0->y();

    double x1 = V1->x();
    double y1 = V1->y();

    double dx01 = x0 - x1;
    double dy01 = y0 - y1;

    double d01 = sqrt(dx01 * dx01 + dy01 * dy01);

    double xs = SymmetryLine->start()->x();
    double ys = SymmetryLine->start()->y();

    double xe = SymmetryLine->end()->x();
    double ye = SymmetryLine->end()->y();

    double dxes = xe - xs;
    double dyes = ye - ys;

    double dse = sqrt(dxes * dxes + dyes * dyes);

    double scale = 1.0 / sqrt(d01 * dse);

    r(EquationIndex) = ((x0 - x1) * (xe - xs) + (y0 - y1) * (ye - ys)) * scale;

    J(EquationIndex, SymmetryLine->start()->X->get_index()) -= dx01 * scale;
    J(EquationIndex, SymmetryLine->start()->Y->get_index()) -= dy01 * scale;

    J(EquationIndex, SymmetryLine->end()->X->get_index()) += dx01 * scale;
    J(EquationIndex, SymmetryLine->end()->Y->get_index()) += dy01 * scale;

    J(EquationIndex, V0->X->get_index()) += dxes * scale;
    J(EquationIndex, V0->Y->get_index()) += dyes * scale;

    J(EquationIndex, V1->X->get_index()) -= dxes * scale;
    J(EquationIndex, V1->Y->get_index()) -= dyes * scale;

    r(EquationIndex + 1) =
            (((x0 - xs) * (y0 - ye) - (x0 - xe) * (y0 - ys)) + ((x1 - xs) * (y1 - ye) - (x1 - xe) * (y1 - ys))) * scale;

    J(EquationIndex + 1, SymmetryLine->start()->X->get_index()) -= (y1 + y0 - 2.0 * ye) * scale;
    J(EquationIndex + 1, SymmetryLine->start()->Y->get_index()) += (x0 + x1 - 2.0 * xe) * scale;

    J(EquationIndex + 1, SymmetryLine->end()->X->get_index()) += (y1 + y0 - 2.0 * ye) * scale;
    J(EquationIndex + 1, SymmetryLine->end()->Y->get_index()) -= (x0 + x1 - 2.0 * xe) * scale;

    J(EquationIndex + 1, V0->X->get_index()) -= dyes * scale;
    J(EquationIndex + 1, V0->Y->get_index()) += dxes * scale;

    J(EquationIndex + 1, V1->X->get_index()) -= dyes * scale;
    J(EquationIndex + 1, V1->Y->get_index()) += dxes * scale;
}