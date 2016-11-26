#ifndef OERSTED_DISTANCE_H
#define OERSTED_DISTANCE_H

#include "Constraint.h"

template<class T>
class Distance : public Constraint {
public:
    Distance(std::shared_ptr<T> e0, std::shared_ptr<T> e1, double d) : Element0(e0), Element1(e1), Dim(d) {};

    size_t set_equation_index(size_t i) override;

    double dim() const { return Dim; };

    void dim(double d) { Dim = d; };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override;

protected:
    std::shared_ptr<T> Element0;
    std::shared_ptr<T> Element1;

    double Dim;
};

#endif //OERSTED_DISTANCE_H
