#ifndef VERTEX_H
#define VERTEX_H

#include "Sketch.h"

class Vertex : public SketchElement {
public:
    Variable *X, *Y;

    // Double Constructors
    Vertex() : X(new Variable(0.0)), Y(new Variable(0.0)) {};

    Vertex(double x, double y) : X(new Variable(x)), Y(new Variable(y)) {};

    Vertex(std::pair<double, double> &xy) : Vertex(xy.first, xy.second) {};

    // Variable Constructors
    Vertex(Variable *xin, Variable *yin) : X(xin), Y(yin) {};

    Vertex(const Vertex &v) {
        X = v.X;
        Y = v.Y;
    };

    Vertex(Variable &x, Variable &y) : X(&x), Y(&y) {};

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

    inline bool operator==(const Vertex &v) { return (v.X == X) && (v.Y == Y); };

    const double x() const { return X->value(); };

    const double y() const { return Y->value(); };

    double hypot() const { return std::hypot(x(), y()); };

    double atan() const { return std::atan2(y(), x()); };

    std::pair<double, double> rotate(const Vertex *origin, const double angle) const;
};

#endif