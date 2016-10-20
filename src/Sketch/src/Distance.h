#ifndef OERSTED_DISTANCE_H
#define OERSTED_DISTANCE_H

template<class T>
class Distance : public Constraint {
public:
    T *Element0;
    T *Element1;

    double Dim;

    // Constructors
    Distance(T &e0, T &e1, double d) : Element0(&e0), Element1(&e1), Dim(d) {};

    // Public Member Functions
    size_t set_equation_index(size_t i) override;

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;
};

#endif //OERSTED_DISTANCE_H
