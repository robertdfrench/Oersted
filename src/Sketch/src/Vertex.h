#ifndef OERSTED_VERTEX_H
#define OERSTED_VERTEX_H

#include "SketchElement.h"
#include "Variable.h"

class Vertex : public SketchElement {
public:
    Vertex() : X(std::make_shared<Variable>(0.0)), Y(std::make_shared<Variable>(0.0)) {};

    Vertex(double x, double y) : X(std::make_shared<Variable>(x)), Y(std::make_shared<Variable>(y)) {};

    Vertex(std::pair<double, double> &xy) : Vertex(xy.first, xy.second) {};

    Vertex(std::shared_ptr<Variable> x, std::shared_ptr<Variable> y) : X(x), Y(y) {};

    Vertex(std::shared_ptr<Vertex> &v) : X(v->X), Y(v->Y) {};

    size_t set_equation_index(size_t i) override { //
        EquationIndex = i;
        return 0;
    };

    size_t x_index() const { return X->get_index();};

    size_t y_index() const { return Y->get_index();};

    bool operator==(Vertex const &v) { return (v.X == X) && (v.Y == Y); };

    double atan() const { return std::atan2(y(), x()); };

    double hypot() const { return std::hypot(x(), y()); };

    double x() const { return X->value(); };

    double y() const { return Y->value(); };

    void register_parameters(Sketch *s) override;

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override {};

    std::pair<double, double> rotate(std::shared_ptr<Vertex> const &origin, const double angle) const;

protected:
    std::shared_ptr<Variable> X;
    std::shared_ptr<Variable> Y;
};

#endif //OERSTED_VERTEX_H