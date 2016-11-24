#ifndef OERSTED_VERTEX_H
#define OERSTED_VERTEX_H

#include "Sketch.h"

class Vertex : public SketchElement {
public:
    std::shared_ptr<Variable> X;
    std::shared_ptr<Variable> Y;

    // Double Constructors
    Vertex() : X(std::make_shared<Variable>(0.0)), Y(std::make_shared<Variable>(0.0)) {};

    Vertex(double x, double y) : X(std::make_shared<Variable>(x)), Y(std::make_shared<Variable>(y)) {};

    Vertex(std::pair<double, double> &xy) : Vertex(xy.first, xy.second) {};

    // Variable Constructors
    Vertex(std::shared_ptr<Variable> x, std::shared_ptr<Variable> y) : X(x), Y(y) {};

    Vertex(std::shared_ptr<Vertex> &v) : X(v->X), Y(v->Y) {};

    //Vertex(Variable &x, Variable &y) : X(&x), Y(&y) {};

    //Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 0;
    };

    void register_parameters(Sketch *s) override {
        s->add_parameter(X);
        s->add_parameter(Y);
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override {};

    bool operator==(Vertex const &v) { return (v.X == X) && (v.Y == Y); };

    const double x() const { return X->value(); };

    const double y() const { return Y->value(); };

    double hypot() const { return std::hypot(x(), y()); };

    double atan() const { return std::atan2(y(), x()); };

    std::pair<double, double> rotate(std::shared_ptr<Vertex> const &origin, const double angle) const;
};

#endif //OERSTED_VERTEX_H