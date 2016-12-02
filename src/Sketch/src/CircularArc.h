#ifndef OERSTED_CIRCULARARC_H
#define OERSTED_CIRCULARARC_H

#include "Curve.h"
#include "Sketch.h"

class CircularArc final : public Curve {
public:
    using Curve::on_manifold;

    using Curve::on_segment;

public:
    CircularArc() : Curve(), Radius(std::make_shared<Variable const>(0.0)) {};

    CircularArc(CircularArc const *c) : Curve(c->Start, c->End, c->ForConstruction), Center(c->Center), Radius(c->Radius) {};

    CircularArc(std::shared_ptr<Vertex const> v0, std::shared_ptr<Vertex const> v1, std::shared_ptr<Vertex const> c, bool fc = false) : Curve(v0, v1, fc), Center(c) {
        double xc = c->x();
        double yc = c->y();
        double x0 = v0->x();
        double y0 = v0->y();
        double x1 = v1->x();
        double y1 = v1->y();

        double r0 = sqrt((xc - x0) * (xc - x0) + (yc - y0) * (yc - y0));
        double r1 = sqrt((xc - x1) * (xc - x1) + (yc - y1) * (yc - y1));

        Radius = std::make_shared<Variable>((r0 + r1) / 2.0);
    };

    CircularArc(std::shared_ptr<Vertex const> v0, std::shared_ptr<Vertex const> v1, std::shared_ptr<Vertex const> c, double r, bool fc = false) : Curve(v0, v1, fc), Center(c), Radius(std::make_shared<Variable const>(r)) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    size_t radius_index() const { return Radius->get_index(); };

    bool is_coincident(std::shared_ptr<Curve const> const &c) const override;

    double a(double s, bool orientation) const override;

    double area() const override;

    double da(double s, bool orientation) const override;

    double length() const override { return radius() * arc_angle(); };

    double radius() const { return Radius->value(); };

    void get_verticies(std::list<std::shared_ptr<Vertex const>> &v, MatchOrientation dir = MatchOrientation::Forward) const override {
        if (dir == MatchOrientation::Forward) {
            v.push_back(Start);
            v.push_back(End);
            v.push_back(Center);
        } else if (dir == MatchOrientation::Reverse) {
            v.push_back(End);
            v.push_back(Start);
            v.push_back(Center);
        }
    };

    void register_parameters(Sketch *s) const override { s->add_parameter(std::const_pointer_cast<Variable>(Radius)); };

    void replace_verticies(std::vector<std::shared_ptr<Vertex const>> const &oldv, std::vector<std::shared_ptr<Vertex const>> const &newv) override;

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override;

    std::shared_ptr<Curve> clone() const override { return std::make_shared<CircularArc>(this); };

    std::shared_ptr<Vertex const> center() const { return Center; };

    MatchOrientation is_identical(std::shared_ptr<Curve const> const &c) const override;

    MatchOrientation is_identical(std::shared_ptr<Curve const> const &c, std::shared_ptr<Vertex const> const &origin, double angle) const override;

    double2 point(double s) const override;

    double2 supremum() const override;

    double2 tangent(double s, bool orientation) const override;

protected:
    std::shared_ptr<Vertex const> Center;
    std::shared_ptr<Variable const> Radius;

    bool on_manifold(double x, double y) const override;

    bool on_segment(double x, double y) const override;

    double s_to_a(double s) const;

    double a_to_s(double a) const;

    double arc_angle() const;

    MatchOrientation is_identical(double r, double xc, double yc, double xs, double ys, double xe, double ye) const;
};

#endif //OERSTED_CIRCULARARC_H
