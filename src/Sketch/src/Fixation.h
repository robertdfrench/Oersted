#ifndef OERSTED_FIXATION_H
#define OERSTED_FIXATION_H

class Fixation : public Constraint {
public:
    Vertex *Point;
    Vertex *Dim;

    // Constructors
    Fixation(Vertex &v);

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

#endif //OERSTED_FIXATION_H