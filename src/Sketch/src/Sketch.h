#ifndef OERSTED_SKETCH_H
#define OERSTED_SKETCH_H

#define SIGN(x) (double)((x > 0.0) - (x < 0.0))

#include <cmath>
#include <vector>
#include <list>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cfloat>
#include "Eigen"

#include <boost/filesystem.hpp>

class Sketch;

// Sketch Parameter
class SketchParameter {
public:
    SketchParameter() : Index(SIZE_MAX) {};

    virtual size_t set_index(size_t i) = 0;

    size_t get_index() const { return Index; };

    virtual void update(Eigen::VectorXd &delta) = 0;

protected:
    size_t Index;
};

class Variable : public SketchParameter {
public:
    Variable() : SketchParameter() {};

    Variable(double v) : SketchParameter(), Value(v) {};

    size_t set_index(size_t i) override {
        Index = i;
        return 1;
    };

    void update(Eigen::VectorXd &delta) override { Value -= delta(Index); };

    const double value() { return Value; };

private:
    double Value;
};

// Sketch Element
class SketchElement {
public:
    virtual size_t set_equation_index(size_t i) = 0;

    size_t get_equation_index() const { return EquationIndex; };

    virtual void register_parameters(Sketch *s) = 0;

    virtual void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) = 0;

protected:
    size_t EquationIndex;
};

class Vertex;

class Curve;

class CircularArc;

class LineSegment;

class Constraint;

class Angle;

template<class T>
class Coincident;

template<class T>
class Distance;

class Fixation;

class Horizontal;

class Length;

class Parallel;

class Perpendicular;

class Radius;

class Rotation;

class Symmetry;

class Tangency;

class Vertical;

class Pattern;

class MirrorCopy;

class RotateCopy;

// Contour
class Contour;

// Utility
struct Branch;

class Star;

class Consellation;

// Sketch
enum class SaveMethod {
    Rasterize
};

class Sketch {
public:
    // Constructors
    Sketch();

    // Public Member Functions
    const Variable *variable(size_t i) const { return Variables[i]; };

    const Vertex *vertex(size_t i) const { return Verticies[i]; };

    const Curve *curve(size_t i) const { return Curves[i]; };

    const Constraint *constraint(size_t i) const { return Constraints[i]; };

    const Contour *contour(size_t i) const { return Contours[i]; };

    const Contour *boundary() const { return Boundary; };

    size_t size() const { return (Verticies.size() + Curves.size() + Constraints.size()); };

    size_t size_verticies() const { return Verticies.size(); };

    size_t size_curves() const { return Curves.size(); };

    size_t size_constraints() const { return Constraints.size(); };

    size_t size_contours() const { return Contours.size(); };

    template<typename T, typename...Args>
    T *new_parameter(Args &&... args);

    template<typename T, typename...Args>
    T &new_element(Args &&... args);

    //#TODO: Add elements by pointer
    void add_element(Vertex &v);

    void add_element(Curve &c);

    void add_element(Constraint &c);

    void add_element(Pattern &p);

    void solve();

    bool build();

    void add_parameter(Variable *v) { add_parameter(Variables, v); };

    // Save Functions
    template<SaveMethod S>
    void save_as(std::string path, std::string file_name) const;

private:
    // SketchParameter
    template<class T>
    void add_parameter(std::vector<T *> &parameters, T *p);

    template<class T>
    void add_element(std::vector<T *> &elements, T &e);

    std::vector<Variable *> Variables;

    std::vector<Vertex *> Verticies;
    std::vector<Curve *> Curves;
    std::vector<Constraint *> Constraints;
    std::vector<Pattern *> Patterns;

    std::vector<Contour *> Contours;

    Contour *Boundary;

    size_t NumVariables;
    size_t NumEquations;
};

// Template Member Functions
// #TODO: Pass elements by pointer
template<typename T>
void Sketch::add_element(std::vector<T *> &elements, T &e) {
    bool unique = true;
    for (size_t i = 0; i < elements.size() && unique; ++i) {
        unique = (elements[i] != &e);
    }

    if (unique) {
        NumEquations += e.set_equation_index(NumEquations);
        elements.push_back(&e);
    }

    e.register_parameters(this);
}

template<typename T>
void Sketch::add_parameter(std::vector<T *> &parameters, T *p) {
    bool unique = true;
    for (size_t i = 0; i < parameters.size() && unique; ++i) {
        unique = (parameters[i] != p);
    }

    if (unique) {
        NumVariables += p->set_index(NumVariables);
        parameters.push_back(p);
    }
}

template<typename T, typename...Args>
T &Sketch::new_element(Args &&... args) {
    T *e = new T(std::forward<Args>(args)...);
    add_element(*e);
    return *e;
}

template<typename T, typename...Args>
T *Sketch::new_parameter(Args &&... args) {
    T *p = new T(std::forward<Args>(args)..., this);
    add_parameter(p);

    return p;
}

#endif //OERSTED_SKETCH_H