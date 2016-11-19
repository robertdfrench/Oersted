#ifndef OERSTED_HORIZONTAL_H
#define OERSTED_HORIZONTAL_H

class Horizontal : public Constraint {
public:
    std::shared_ptr<LineSegment> Line;

    Horizontal(std::shared_ptr<LineSegment> l) : Line(l) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

#endif //OERSTED_HORIZONTAL_H
