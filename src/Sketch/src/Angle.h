#ifndef OERSTED_ANGLE_H
#define OERSTED_ANGLE_H

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


#endif //OERSTED_ANGLE_H
