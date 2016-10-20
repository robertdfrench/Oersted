#ifndef OERSTED_ROTATION_H
#define OERSTED_ROTATION_H

class Rotation : public Constraint {
public:
    Vertex *V0;
    Vertex *V1;
    Vertex *Origin;
    double Angle;

    Rotation(Vertex &v0, Vertex &v1, Vertex &origin, double a) : V0(&v0), V1(&v1), Origin(&origin), Angle(a) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r);
};

#endif //OERSTED_ROTATION_H
