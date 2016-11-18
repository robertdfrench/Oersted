#ifndef OERSTED_ROTATION_H
#define OERSTED_ROTATION_H

class Rotation : public Constraint {
public:
    std::shared_ptr<Vertex> V0;
    std::shared_ptr<Vertex> V1;
    std::shared_ptr<Vertex> Origin;
    double Angle;

    Rotation(std::shared_ptr<Vertex> v0, std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> origin, double a) : V0(v0), V1(v1), Origin(origin), Angle(a) {};

    size_t set_equation_index(size_t i) override {
        EquationIndex = i;
        return 2;
    };

    void update(Eigen::MatrixXd &J, Eigen::VectorXd &r);
};

#endif //OERSTED_ROTATION_H
