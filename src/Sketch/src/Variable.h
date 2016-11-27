#ifndef OERSTED_VARIABLE_H
#define OERSTED_VARIABLE_H

#include <cfloat>

#include "Eigen"

class Variable {
public:
    Variable() : Value(DBL_MAX), Index(SIZE_MAX) {};

    Variable(double v) : Value(v), Index(SIZE_MAX) {};

    size_t get_index() const { return Index; };

    size_t set_index(size_t i) {
        Index = i;
        return 1;
    };

    double value() const { return Value; };

    void update(Eigen::VectorXd &delta) { Value -= delta(Index); };

protected:
    double Value;

    size_t Index;
};

#endif //OERSTED_VARIABLE_H
