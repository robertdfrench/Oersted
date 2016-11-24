#ifndef OERSTED_CIRCULARARC_H
#define OERSTED_CIRCULARARC_H

#include "Curve.h"

class CircularArc final : public Curve {
public:
    friend class Coincident<CircularArc>;

    friend class Distance<CircularArc>;

    friend class Radius;

    friend class Tangency;

    CircularArc() : Curve(), Radius(std::make_shared<Variable>(0.0)) {};

    CircularArc(CircularArc const *c) : Curve(c->Start, c->End, c->ForConstruction), Center(c->Center), Radius(c->Radius) {};

    CircularArc(std::shared_ptr<Vertex> v0, std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> c, bool fc = false) : Curve(v0, v1, fc), Center(c) {};

    CircularArc(std::shared_ptr<Vertex> v0, std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> c, double r, bool fc = false) : Curve(v0, v1, fc), Center(c), Radius(std::make_shared<Variable>(r)) {};

    CircularArc(std::shared_ptr<Vertex> v0, std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> c, std::shared_ptr<Variable> r, Sketch &s, bool fc = false) : Curve(v0, v1, fc), Center(c), Radius(r) {};

    std::shared_ptr<Vertex> center() const { return Center; };

    double radius() const { return Radius->value(); };

    void get_verticies(std::list<std::shared_ptr<Vertex>> &v, Direction dir = Direction::Forward) const override {
        if (dir == Direction::Forward) {
            v.push_back(Start);
            v.push_back(End);
            v.push_back(Center);
        } else if (dir == Direction::Reverse) {
            v.push_back(End);
            v.push_back(Start);
            v.push_back(Center);
        }
    };

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void register_parameters(Sketch *s) override { s->add_parameter(Radius); };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;

    sPoint point(double s) const override;

    Vertex tangent(double s, bool orientation) const override;

    double length() const override { throw; };

    double area() const override;

    double a(double s, bool orientation) const override;

    double da(double s, bool orientation) const override;

    std::pair<double, double> supremum() const override;

    using Curve::on_manifold;

    using Curve::on_segment;

    Direction is_identical(std::shared_ptr<Curve> const &c) const override;

    Direction is_identical(std::shared_ptr<Curve> const &c, std::shared_ptr<Vertex> const &origin, double const angle) const override;

    bool is_coincident(std::shared_ptr<Curve> const &c) const override;

    std::shared_ptr<Curve> clone() const override { return std::make_shared<CircularArc>(this); };

    void replace_verticies(std::vector<std::shared_ptr<Vertex>> oldv, std::vector<std::shared_ptr<Vertex>> newv) override;

protected:
    std::shared_ptr<Vertex> Center;
    std::shared_ptr<Variable> Radius;

    bool on_manifold(const double x, const double y) const override;

    bool on_segment(const double x, const double y) const override;

    Direction is_identical(const double r, const double xc, const double yc, const double xs, const double ys, const double xe, const double ye) const;

private:
    double s_to_a(double s) const;

    double a_to_s(double a) const;

    double arc_angle() const;
};

#endif //OERSTED_CIRCULARARC_H
