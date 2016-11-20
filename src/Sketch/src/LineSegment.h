#ifndef OERSTED_LINESEGMENT_H
#define OERSTED_LINESEGMENT_H

#include "Curve.h"

class LineSegment final : public Curve { // TODO: Use interface instead of friend class
public:
    friend class Angle;

    friend class Coincident<LineSegment>;

    friend class Distance<LineSegment>;

    friend class Horizontal;

    friend class Length;

    friend class Tangency;

    friend class Vertical;

    //Constructors
    LineSegment() : Curve() {};

    LineSegment(LineSegment const *l) : Curve(l->Start, l->End, l->ForConstruction) {};

    LineSegment(std::shared_ptr<Vertex> v0, std::shared_ptr<Vertex> v1, bool fc = false) : Curve(v0, v1, fc) {};

    // Virtual Function Implementation
    void get_verticies(std::list<std::shared_ptr<Vertex>> &v) const override {
        v.push_back(Start);
        v.push_back(End);
    };

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 0;
    };

    void register_parameters(Sketch *s) override {};

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override {};

    // Calculation
    sPoint point(double s) const override;

    Vertex tangent(double s, bool orientation) const override;

    double length() const override;

    double area() const override { return 0.0; };

    double a(double s, bool orientation) const override;

    double da(double s, bool orientation) const override { return 0.0; };

    std::pair<double, double> supremum() const override;

    // Curve-Vertex Comparison
    using Curve::on_manifold;

    using Curve::on_segment;

    // Curve-Curve Comparison
    bool is_identical(std::shared_ptr<Curve> const &c) const override;

    bool is_identical(std::shared_ptr<Curve> const &c, std::shared_ptr<Vertex> const &origin, double const angle) const override;

    bool is_coincident(std::shared_ptr<Curve> const &c) const override;

    // Modification
    std::shared_ptr<Curve> clone() const override { return std::make_shared<LineSegment>(this); };

    void replace_verticies(std::vector<std::shared_ptr<Vertex>> oldv, std::vector<std::shared_ptr<Vertex>> newv) override;

protected:
    bool on_manifold(const double x, const double y) const override;

    bool on_segment(const double x, const double y) const override;

    bool is_identical(const double x0, const double y0, const double x1, const double y1) const;
};

#endif //OERSTED_LINESEGMENT_H
