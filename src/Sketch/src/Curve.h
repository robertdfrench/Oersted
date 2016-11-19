#ifndef OERSTED_CURVE_H
#define OERSTED_CURVE_H

#include "Sketch.h"

class Curve : public SketchElement {
public:
    friend class MirrorCopy;

    // Constructors
    Curve() : Start(nullptr), End(nullptr) {};

    Curve(std::shared_ptr<Vertex> v0, std::shared_ptr<Vertex>v1, bool fc = false) : Start(v0), End(v1), ForConstruction(fc) {};

    //Curve(Vertex *v0, Vertex *v1, bool fc = false) : Start(v0), End(v1), ForConstruction(fc) {};

    // Properties
    bool ForConstruction = false;

    // Accessors
    std::shared_ptr<Vertex> start() const { return Start; };

    std::shared_ptr<Vertex> end() const { return End; };

    virtual void get_verticies(std::list<std::shared_ptr<Vertex>> &v) const = 0;

    // Calculation
    // #TODO: Add const to double_t and bool arguments where appropriate
    virtual sPoint point(double s) const = 0;

    virtual Vertex tangent(double s, bool orientation) const = 0;

    virtual double length() const = 0; // length of segment between end points
    virtual double area() const = 0;    // area of segment between end points

    virtual double a(double s, bool orientation) const = 0;     // tangent angle
    virtual double da(double s, bool orientation) const = 0;    // derivative of tangent angle per unit length

    virtual std::pair<double, double> supremum() const = 0;     // maximum length of vector between origin and point on curve

    // Curve-Vertex Comparison
    virtual bool on_manifold(std::shared_ptr<Vertex> v) const final; // true if vertex is on manifold defined by curve
    virtual bool on_manifold(std::shared_ptr<Vertex> v, std::shared_ptr<Vertex> origin, const double angle) const final;

    virtual bool on_segment(std::shared_ptr<Vertex> v) const final; // true if vertex is on curve segment
    virtual bool on_segment(std::shared_ptr<Vertex> v, std::shared_ptr<Vertex> origin, const double angle) const final;

    // Curve-Curve Comparison
    virtual bool is_identical(std::shared_ptr<Curve> c) const = 0; // true if (input curve) XOR (object curve) is a set with measure < tol
    virtual bool is_identical(std::shared_ptr<Curve> c, std::shared_ptr<Vertex> origin, const double angle) const = 0;

    // #TODO: virtual bool is_overlapping(const Curve* c) const = 0; // true if (input curve) AND (object curve) is a set with measure > tol
    // #TODO: virtual bool is_overlapping(const Curve* c, const Vertex* origin, const double_t angle) const = 0;

    virtual bool is_coincident(std::shared_ptr<Curve> c) const = 0; // true if (input curve) AND (object curve + parametric extension) is a set with measure > tol
    // #TODO: virtual bool is_coincident(const Curve* c, const Vertex* origin, const double_t angle) const = 0;

    // Modification
    void reverse() { std::swap(Start, End); };

    virtual std::shared_ptr<Curve> clone() const = 0;

    virtual void replace_verticies(std::vector<std::shared_ptr<Vertex>> oldv, std::vector<std::shared_ptr<Vertex>> newv) = 0;

protected:
    std::shared_ptr<Vertex> Start;
    std::shared_ptr<Vertex> End;

    virtual bool on_manifold(const double x, const double y) const = 0;

    virtual bool on_segment(const double x, const double y) const = 0;
};

#include "LineSegment.h"
#include "CircularArc.h"

#endif //OERSTED_CURVE_H