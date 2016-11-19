#ifndef OERSTED_SYMMETRY_H
#define OERSTED_SYMMETRY_H

class Symmetry : public Constraint {
public:
    std::shared_ptr<Vertex> V0;
    std::shared_ptr<Vertex> V1;
    std::shared_ptr<LineSegment> SymmetryLine;

    Symmetry(std::shared_ptr<Vertex> v0, std::shared_ptr<Vertex> v1, std::shared_ptr<LineSegment> line) : V0(v0), V1(v1), SymmetryLine(line) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

#endif //OERSTED_SYMMETRY_H
