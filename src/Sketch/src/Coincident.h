#ifndef OERSTED_COINCIDENT_H
#define OERSTED_COINCIDENT_H

template<class T>
class Coincident : public Constraint {
public:
    Vertex *Point;
    T *Element;

    // Constructors
    Coincident(Vertex &p, T &e) : Point(&p), Element(&e) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

#endif //OERSTED_COINCIDENT_H
