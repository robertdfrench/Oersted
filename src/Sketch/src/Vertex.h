#ifndef OERSTED_VERTEX_H
#define OERSTED_VERTEX_H

#include "SketchElement.h"
#include "Variable.h"
#include "doublen.h"

class Vertex : public SketchElement {
public:
    Vertex() : X(std::make_shared<Variable const>()), Y(std::make_shared<Variable const>()) {};

    Vertex(double x, double y) : X(std::make_shared<Variable const>(x)), Y(std::make_shared<Variable const>(y)) {};

    Vertex(std::shared_ptr<Variable const> x, std::shared_ptr<Variable const> y) : X(x), Y(y) {};

    Vertex(std::shared_ptr<Vertex const> const &v) : X(v->X), Y(v->Y) {};

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

    void register_parameters(Sketch *s) const override;

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override {};

    double2 rotate(std::shared_ptr<Vertex const> const &origin, double angle) const;

protected:
    std::shared_ptr<Variable const> X;

    std::shared_ptr<Variable const> Y;
};

#endif //OERSTED_VERTEX_H