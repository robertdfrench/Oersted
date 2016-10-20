#ifndef OERSTED_CURVE_H
#define OERSTED_CURVE_H

#include "Sketch.h"

class Curve : public SketchElement {
public:
    friend class MirrorCopy;

    // Constructors
    Curve() : Start(nullptr), End(nullptr) {};

    Curve(Vertex &v0, Vertex &v1) : Start(&v0), End(&v1) {};

    Curve(Vertex *v0, Vertex *v1) : Start(v0), End(v1) {};

    // Properties
    bool ForConstruction = false;

    // Accessors
    const Vertex *start() const { return Start; };

    const Vertex *end() const { return End; };

    virtual void get_verticies(std::list<Vertex *> &v) const = 0;

    // Calculation
    // #TODO: Add const to double_t and bool arguments where appropriate
    virtual Vertex point(double s) const = 0;

    virtual Vertex tangent(double s, bool orientation) const = 0;

    virtual double length() const = 0; // length of segment between end points
    virtual double area() const = 0;    // area of segment between end points

    virtual double a(double s, bool orientation) const = 0;        // tangent angle
    virtual double da(double s, bool orientation) const = 0;    // derivative of tangent angle per unit length

    virtual double supremum() const = 0;    // maximum length of vector between origin and point on curve

    // Curve-Vertex Comparison
    virtual bool on_manifold(const Vertex *v) const final; // true if vertex is on manifold defined by curve
    virtual bool on_manifold(const Vertex *v, const Vertex *origin, const double angle) const final;

    virtual bool on_segment(const Vertex *v) const final; // true if vertex is on curve segment
    virtual bool on_segment(const Vertex *v, const Vertex *origin, const double angle) const final;

    // Curve-Curve Comparison
    virtual bool is_identical(const Curve *c) const = 0; // true if (input curve) XOR (object curve) is a set with measure < tol
    virtual bool is_identical(const Curve *c, const Vertex *origin, const double angle) const = 0;

    // #TODO: virtual bool is_overlapping(const Curve* c) const = 0; // true if (input curve) AND (object curve) is a set with measure > tol
    // #TODO: virtual bool is_overlapping(const Curve* c, const Vertex* origin, const double_t angle) const = 0;

    virtual bool is_coincident(const Curve *c) const = 0; // true if (input curve) AND (object curve + parametric extension) is a set with measure > tol
    // #TODO: virtual bool is_coincident(const Curve* c, const Vertex* origin, const double_t angle) const = 0;

    // Modification
    Curve *split(Vertex *v, double s);

    void reverse() { std::swap(Start, End); };

    virtual Curve *clone() const = 0;

    virtual void replace_verticies(std::vector<Vertex *> oldv, std::vector<Vertex *> newv) = 0;

protected:
    Vertex *Start, *End;

    virtual bool on_manifold(const double x, const double y) const = 0;

    virtual bool on_segment(const double x, const double y) const = 0;
};

#include "LineSegment.h"
#include "CircularArc.h"

#endif //OERSTED_CURVE_H