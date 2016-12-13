#ifndef OERSTED_COINCIDENT_H
#define OERSTED_COINCIDENT_H

#include "Constraint.h"

class Vertex;

template<class T>
class Coincident : public Constraint {
public:
    Coincident(std::shared_ptr<Vertex const> p, std::shared_ptr<T const> e) : Point(p), Element(e) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 1;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) const override;

protected:
    std::shared_ptr<Vertex const> Point;

    std::shared_ptr<T const> Element;
};

#endif //OERSTED_COINCIDENT_H
