#ifndef OERSTED_SYMMETRY_H
#define OERSTED_SYMMETRY_H

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

#endif //OERSTED_SYMMETRY_H
