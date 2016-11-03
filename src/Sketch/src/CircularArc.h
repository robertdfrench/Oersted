#ifndef OERSTED_CIRCULARARC_H
#define OERSTED_CIRCULARARC_H

#include "Curve.h"

class CircularArc final : public Curve {
public:
    friend class Coincident<CircularArc>;

    friend class Distance<CircularArc>;

    friend class Radius;

    friend class Tangency;

    // Constructors
    CircularArc() : Curve(), Radius(new Variable(0.0)) {};

    CircularArc(const CircularArc *c) : Curve(c->Start, c->End, c->ForConstruction), Center(c->Center), Radius(c->Radius) {};

    CircularArc(Vertex &v0, Vertex &v1, Vertex &c, bool fc = false) : Curve(v0, v1, fc), Center(&c) {};

    CircularArc(Vertex &v0, Vertex &v1, Vertex &c, double r, bool fc = false) : Curve(v0, v1, fc), Center(&c), Radius(new Variable{r}) {};

    CircularArc(Vertex &v0, Vertex &v1, Vertex &c, Variable &r, Sketch &s, bool fc = false) : Curve(v0, v1, fc), Center(&c), Radius(&r) {};

    // Accessors
    const Vertex *center() const { return Center; };

    double radius() const { return Radius->value(); };

    // Virtual Function Implementation
    void get_verticies(std::list<Vertex *> &v) const override {
        v.push_back(Start);
        v.push_back(End);
        v.push_back(Center);
    };

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void register_parameters(Sketch *s) override { s->add_parameter(Radius); };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;

    // Calculation
    Vertex point(double s) const override;

    Vertex tangent(double s, bool orientation) const override;

    double length() const override { throw; };

    double area() const override;

    double a(double s, bool orientation) const override;

    double da(double s, bool orientation) const override;

    std::pair<double, double> supremum() const override;

    // Curve-Vertex Comparison
    using Curve::on_manifold;

    using Curve::on_segment;

    // Curve-Curve Comparison
    bool is_identical(const Curve *c) const override;

    bool is_identical(const Curve *c, const Vertex *origin, const double angle) const override;

    bool is_coincident(const Curve *c) const override;

    // Modification
    Curve *clone() const override { return new CircularArc(this); };

    void replace_verticies(std::vector<Vertex *> oldv, std::vector<Vertex *> newv) override;

protected:
    Vertex *Center;
    Variable *Radius;

    bool on_manifold(const double x, const double y) const override;

    bool on_segment(const double x, const double y) const override;

    bool is_identical(const double r, const double xc, const double yc, const double xs, const double ys, const double xe, const double ye) const;

private:
    double s_to_a(double s) const;

    double a_to_s(double a) const;

    double arc_angle() const;
};

#endif //OERSTED_CIRCULARARC_H
