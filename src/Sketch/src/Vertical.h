#ifndef OERSTED_VERTICAL_H
#define OERSTED_VERTICAL_H

class Vertical : public Constraint {
public:
    std::shared_ptr<LineSegment> Line;

    // Constructors
    Vertical(std::shared_ptr<LineSegment> l) : Line(l) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

#endif //OERSTED_VERTICAL_H
