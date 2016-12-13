#ifndef OERSTED_LINESEGMENT_H
#define OERSTED_LINESEGMENT_H

#include "Curve.h"
#include "Sketch.h"

class LineSegment final : public Curve {
public:
    using Curve::on_manifold;
    using Curve::on_segment;

public:
    LineSegment() : Curve() {};

    LineSegment(LineSegment const *l) : Curve(l->Start, l->End, l->ForConstruction) {};

    LineSegment(std::shared_ptr<Vertex const> v0, std::shared_ptr<Vertex const> v1, bool fc = false) : Curve(v0, v1, fc) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 0;
    };

    bool is_coincident(std::shared_ptr<Curve const> const &c) const override;

    double length() const override;

    double area() const override { return 0.0; };

    double a(double s, bool orientation) const override;

    double da(double s, bool orientation) const override { return 0.0; };

    void get_verticies(std::list<std::shared_ptr<Vertex const>> &v, MatchOrientation dir = MatchOrientation::Forward) const override {
        if (dir == MatchOrientation::Forward) {
            v.push_back(Start);
            v.push_back(End);
        } else if (dir == MatchOrientation::Reverse) {
            v.push_back(End);
            v.push_back(Start);
        }
    };

    void register_parameters(Sketch *s) const override {};

    void replace_verticies(std::vector<std::shared_ptr<Vertex const>> const &oldv, std::vector<std::shared_ptr<Vertex const>> const &newv) override;

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override {};

    std::shared_ptr<Curve> clone() const override { return std::make_shared<LineSegment>(this); };

    MatchOrientation is_identical(std::shared_ptr<Curve const> const &c) const override;

    MatchOrientation is_identical(std::shared_ptr<Curve const> const &c, std::shared_ptr<Vertex const> const &origin, double angle) const override;

    double2 point(double s) const override;

    double2 supremum() const override;

    double2 tangent(double s, bool orientation) const override;

protected:
    bool on_manifold(double x, double y) const override;

    bool on_segment(double x, double y) const override;

    MatchOrientation is_identical(double x0, double y0, double x1, double y1) const;
};

#endif //OERSTED_LINESEGMENT_H
