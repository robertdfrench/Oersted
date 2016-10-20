#ifndef OERSTED_LENGTH_H
#define OERSTED_LENGTH_H

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

#endif //OERSTED_LENGTH_H