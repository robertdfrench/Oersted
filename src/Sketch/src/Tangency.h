#ifndef OERSTED_TANGENCY_H
#define OERSTED_TANGENCY_H

class Tangency : public Constraint {
public:
    std::shared_ptr<CircularArc> Arc;
    std::shared_ptr<LineSegment>Line;

    // Constructors
    Tangency(std::shared_ptr<CircularArc> ca, std::shared_ptr<LineSegment> ls) : Arc(ca), Line(ls) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

#endif //OERSTED_TANGENCY_H
