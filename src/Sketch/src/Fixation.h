#ifndef OERSTED_FIXATION_H
#define OERSTED_FIXATION_H

class Fixation : public Constraint {
public:
    std::shared_ptr<Vertex> Point;
    sPoint Dim;

    Fixation(std::shared_ptr<Vertex> v);

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

#endif //OERSTED_FIXATION_H