#ifndef OERSTED_RADIUS_H
#define OERSTED_RADIUS_H

class Radius : public Constraint {
public:
    std::shared_ptr<CircularArc> Arc;

    double Dim;

    Radius(std::shared_ptr<CircularArc> c, double r) : Arc(c), Dim(r) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

#endif //OERSTED_RADIUS_H
