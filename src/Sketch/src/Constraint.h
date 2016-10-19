#ifndef OERSTED_CONSTRAINT_H
#define OERSTED_CONSTRAINT_H

#include "Sketch.h"

class Constraint : public SketchElement {
public:
    void register_parameters(Sketch *s) override {};
};

class Angle : public Constraint {
public:
    LineSegment *Line0, *Line1;

    double Dim;

    // Constructors
    Angle(LineSegment &l0, LineSegment &l1, double angle) : Line0(&l0), Line1(&l1), Dim(angle) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

template<class T>
class Coincident : public Constraint {
public:
    Vertex *Point;
    T *Element;

    // Constructors
    Coincident(Vertex &p, T &e) : Point(&p), Element(&e) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

template<class T>
class Distance : public Constraint {
public:
    T *Element0;
    T *Element1;

    double Dim;

    // Constructors
    Distance(T &e0, T &e1, double d) : Element0(&e0), Element1(&e1), Dim(d) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override;

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

class Fixation : public Constraint {
public:
    Vertex *Point;
    Vertex *Dim;

    // Constructors
    Fixation(Vertex &v);

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

class Horizontal : public Constraint {
public:
    LineSegment *Line;

    // Constructors
    Horizontal(LineSegment &l) : Line(&l) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

class Length : public Constraint {
public:
    LineSegment *Line;

    double Dim;

    // Constructors
    Length(LineSegment &c, double length) : Line(&c), Dim(length) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

class Radius : public Constraint {
public:
    CircularArc *Arc;

    double Dim;

    // Constructors
    Radius(CircularArc &c, double r) : Arc(&c), Dim(r) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

class Rotation : public Constraint {
public:
    Vertex *V0;
    Vertex *V1;
    Vertex *Origin;
    double Angle;

    Rotation(Vertex &v0, Vertex &v1, Vertex &origin, double a) : V0(&v0), V1(&v1), Origin(&origin), Angle(a) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r);
};

class Symmetry : public Constraint {
public:
    Vertex *V0;
    Vertex *V1;
    LineSegment *SymmetryLine;

    Symmetry(Vertex &v0, Vertex &v1, LineSegment &line) : V0(&v0), V1(&v1), SymmetryLine(&line) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

class Tangency : public Constraint {
public:
    CircularArc *Arc;
    LineSegment *Line;

    // Constructors
    Tangency(CircularArc &ca, LineSegment &ls) : Arc(&ca), Line(&ls) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

class Vertical : public Constraint {
public:
    LineSegment *Line;

    // Constructors
    Vertical(LineSegment &l) : Line(&l) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

#endif //OERSTED_CONSTRAINT_H