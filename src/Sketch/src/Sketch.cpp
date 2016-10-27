#include "Sketch.hpp"

Sketch::Sketch() {
    NumEquations = 0;
    NumVariables = 0;
    Boundary = new Contour();
}

void Sketch::solve() {
    // #TODO: Tolerance based iteration convergence monitoring
    // #TODO: Use sparse matrix representation for Jacobian

    Eigen::VectorXd delta = Eigen::VectorXd::Zero(NumVariables);
    Eigen::VectorXd r = Eigen::VectorXd::Zero(NumEquations);
    Eigen::VectorXd dr = Eigen::VectorXd::Zero(NumEquations);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(NumEquations, NumVariables);
    Eigen::MatrixXd JJT = Eigen::MatrixXd::Zero(NumEquations, NumEquations);
    Eigen::FullPivLU<Eigen::MatrixXd> LU;

    for (size_t j = 0; j < 32; ++j) {
        // #TODO: Use QR or LDLT decomposition instead
        // #TODO: Check for over/underdetermined system
        // #TODO: Try extracting sub matrix and solving
        // #TODO: Extract newton_update over Verticies, Curves, and Constraints into a private method

        J.setZero();

        for (size_t i = 0; i < Verticies.size(); ++i) {
            Verticies[i]->update(J, r);
        }

        for (size_t i = 0; i < Curves.size(); ++i) {
            Curves[i]->update(J, r);
        }

        for (size_t i = 0; i < Constraints.size(); ++i) {
            Constraints[i]->update(J, r);
        }

        JJT = J * J.transpose();
        LU = JJT.fullPivLu();
        delta = J.transpose() * LU.solve(r);

        for (size_t i = 0; i < Variables.size(); ++i) {
            Variables[i]->update(delta);
        }
    }
}

bool Sketch::build() {
    Constellation c = Constellation(this);


    bool success = c.boundary(Boundary);

    Contours.resize(0);
    success = c.contours(Contours) && success;

    return success;
}

void Sketch::add_element(Vertex &v) {
    add_element(Verticies, v);
}

void Sketch::add_element(Curve &c) {
    add_element(Curves, c);
}

void Sketch::add_element(Constraint &c) {
    add_element(Constraints, c);
}

void Sketch::add_element(Pattern &p) {
    add_element(Patterns, p);
    p.register_elements(this);
}

template<>
void Sketch::save_as<SaveMethod::Rasterize>(std::string path, std::string file_name) const {
    /*
        This is a stub for visualization
    */

    if(!boost::filesystem::exists(path)) {
        boost::filesystem::create_directories(path);
    }

    std::fstream fs;
    fs.open(path + file_name, std::fstream::out);

    for (size_t i = 0; i < Curves.size(); ++i) {
        if (!Curves[i]->ForConstruction) {
            for (size_t j = 0; j <= 10; ++j) {
                Vertex v = Curves[i]->point(j / 10.0);
                fs << v.x() << "," << v.y() << "\n";
            }
            fs << "NaN" << "," << "NaN" << "\n";
        }
    }

    fs.close();
}