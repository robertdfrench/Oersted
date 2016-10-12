#ifndef LINESEGMENT_H
#define LINESEGMENT_H

#include "Curve.h"

class LineSegment final : public Curve {
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

    LineSegment(const LineSegment *l) : Curve(l->Start, l->End) {};

    LineSegment(Vertex &v0, Vertex &v1) : Curve(v0, v1) {};

    // Virtual Function Implementation
    void get_verticies(std::list<Vertex *> &v) const override {
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
    Vertex point(double s) const override;

    Vertex tangent(double s, bool orientation) const override;

    double length() const override;

    double area() const override { return 0.0; };

    double a(double s, bool orientation) const override;

    double da(double s, bool orientation) const override { return 0.0; };

    double supremum() const override;

    // Curve-Vertex Comparison
    using Curve::on_manifold;

    using Curve::on_segment;

    // Curve-Curve Comparison
    bool is_identical(const Curve *c) const override;

    bool is_identical(const Curve *c, const Vertex *origin, const double angle) const override;

    bool is_coincident(const Curve *c) const override;

    // Modification
    Curve *clone() const override { return new LineSegment(this); };

    void replace_verticies(std::vector<Vertex *> oldv, std::vector<Vertex *> newv) override;

protected:
    bool on_manifold(const double x, const double y) const override;

    bool on_segment(const double x, const double y) const override;

    bool is_identical(const double x0, const double y0, const double x1, const double y1) const;
};

#endif
