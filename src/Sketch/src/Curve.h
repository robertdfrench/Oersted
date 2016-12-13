#ifndef OERSTED_CURVE_H
#define OERSTED_CURVE_H

#include "SketchElement.h"
#include "Vertex.h"

class double2;

enum class MatchOrientation : signed int {
    None = 0,
    Forward = 1,
    Reverse = -1
};

class Curve : public SketchElement {
public:
    Curve() : Start(std::make_shared<Vertex const>()), End(std::make_shared<Vertex const>()), ForConstruction(false) {};

    Curve(std::shared_ptr<Vertex const> v0, std::shared_ptr<Vertex const> v1, bool fc = false) : Start(v0), End(v1), ForConstruction(fc) {};

    bool for_construction() const { return ForConstruction; };

    void for_construction(bool fc) { ForConstruction = fc; };

    void reverse() { std::swap(Start, End); };

    std::shared_ptr<Vertex const> start() const { return Start; };

    std::shared_ptr<Vertex const> end() const { return End; };

public: // interface
    virtual bool is_coincident(std::shared_ptr<Curve const> const &c) const = 0; // true if (input curve) AND (object curve + parametric extension) is a set with measure > tol

    // virtual bool is_coincident(const Curve* c, const Vertex* origin, const double_t angle) const = 0; // TODO

    // virtual bool is_overlapping(const Curve* c) const = 0; // true if (input curve) AND (object curve) is a set with measure > tol, // TODO

    // virtual bool is_overlapping(const Curve* c, const Vertex* origin, const double_t angle) const = 0; // TODO

    virtual bool on_manifold(std::shared_ptr<Vertex const> const &v) const final; // true if vertex is on manifold defined by curve

    virtual bool on_manifold(std::shared_ptr<Vertex const> const &v, std::shared_ptr<Vertex const> const &origin, double angle) const final;

    virtual bool on_segment(std::shared_ptr<Vertex const> const &v) const final; // true if vertex is on curve segment

    virtual bool on_segment(std::shared_ptr<Vertex const> const &v, std::shared_ptr<Vertex const> const &origin, double angle) const final;

    virtual double a(double s, bool orientation) const = 0; // tangent angle

    virtual double area() const = 0; // area of segment between end points

    virtual double da(double s, bool orientation) const = 0; // derivative of tangent angle per unit length

    virtual double length() const = 0; // length of segment between end points

    virtual void get_verticies(std::list<std::shared_ptr<Vertex const>> &v, MatchOrientation dir = MatchOrientation::Forward) const = 0;

    virtual void replace_verticies(std::vector<std::shared_ptr<Vertex const>> const &oldv, std::vector<std::shared_ptr<Vertex const>> const &newv) = 0;

    virtual std::shared_ptr<Curve> clone() const = 0;

    virtual MatchOrientation is_identical(std::shared_ptr<Curve const> const &c) const = 0; // true if (input curve) XOR (object curve) is a set with measure < tol

    virtual MatchOrientation is_identical(std::shared_ptr<Curve const> const &c, std::shared_ptr<Vertex const> const &origin, double const angle) const = 0;

    virtual double2 point(double s) const = 0;

    virtual double2 supremum() const = 0; // maximum length of vector between origin and point on curve

    virtual double2 tangent(double s, bool orientation) const = 0;

protected:
    std::shared_ptr<Vertex const> Start;

    std::shared_ptr<Vertex const> End;

    bool ForConstruction;

    virtual bool on_manifold(double x, double y) const = 0;

    virtual bool on_segment(double x, double y) const = 0;
};

#endif //OERSTED_CURVE_H