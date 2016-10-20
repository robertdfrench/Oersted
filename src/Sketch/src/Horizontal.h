#ifndef OERSTED_HORIZONTAL_H
#define OERSTED_HORIZONTAL_H

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

#endif //OERSTED_HORIZONTAL_H
