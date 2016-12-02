#include <random>

#include "Sketch.h"
#include "Vertex.h"
#include "Curve.h"
#include "Constraint.h"
#include "Pattern.h"
#include "Constellation.h"
#include "Contour.h"
#include "Star.h"
#include "Branch.h"
#include "Variable.h"
#include "doublen.h"

#include <iostream>

double Sketch::solve() {
    // TODO: Time to finish the bells and whistles here
    // TODO: Use sparse matrix representation for Jacobian, (this is lower priority since typical uses will not produce very large matrices)
    // TODO: Add residual instrumentation to tests

    if (NumVariables < NumEquations) {
        std::cerr << "Variables = " << NumVariables << ", Equations = " << NumEquations << ". Number of variables is less than the number of equations. The sketch may be over-constrained." << std::endl;
    }

    Eigen::VectorXd delta = Eigen::VectorXd::Zero(NumVariables);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(NumVariables);
    Eigen::VectorXd r = Eigen::VectorXd::Zero(NumEquations);
    Eigen::VectorXd dr = Eigen::VectorXd::Zero(NumEquations);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(NumEquations, NumVariables);
    //Eigen::MatrixXd JJT = Eigen::MatrixXd::Zero(NumEquations, NumEquations);
    //Eigen::FullPivHouseholderQR<Eigen::MatrixXd> QR;
    //Eigen::FullPivLU<Eigen::MatrixXd> LU;

    // Calculate scale parameter for tolerances
    double scale = 0.0;
    if (Curves.size() > 0) {
        for (size_t i = 0;i!= Curves.size();++i){
            scale = std::max(scale, Curves[i]->length());
        }
    } else if (Verticies.size() > 0) {
        double vmax = DBL_MIN;
        double vmin = DBL_MAX;
        for (size_t i = 0; i != Verticies.size(); ++i) {
            vmax = std::max(vmax, Verticies[i]->x());
            vmin = std::min(vmin, Verticies[i]->x());

            vmax = std::max(vmax, Verticies[i]->y());
            vmin = std::max(vmin, Verticies[i]->y());
        }

        if (vmax != vmin) {
            scale = vmax - vmin;
        } else {
            scale = 1.0;
        }
    }

    // Initial residual calculation
    for (size_t i = 0; i != Verticies.size(); ++i) {
        Verticies[i]->update(J, r);
    }

    for (size_t i = 0; i != Curves.size(); ++i) {
        Curves[i]->update(J, r);
    }

    for (size_t i = 0; i != Constraints.size(); ++i) {
        Constraints[i]->update(J, r);
    }

    // Setup iteration
    double residual_norm{r.norm() / scale};
    double residual_tol = std::sqrt(NumEquations + 1.0) * DBL_EPSILON;
    int armijo_int{0};
    long max_rank(std::min(NumEquations, NumVariables));
    for (size_t j = 0; residual_norm > residual_tol && j != (8 * sizeof(double)); ++j) {
        // TODO: Try extracting sub matrix and solving (Certain constraints are known to eliminate variables (Radius, Fixation), if the matrix can be permuted so that a submatrix on the block diagonal is lower-triangular, then a subsystem of equations can be solved independently of the remaining equations. Moreover, if the jacobian can be permuted into a block lower triangular form, subsets of equations can be solved blockwise)
        // TODO: Extract newton_update over Verticies, Curves, and Constraints into a private method
        // TODO: Extract armijo_update into a private method
        // TODO: Factor out rank deficiency perturbation to be independent of QR/LU/LDL^T choice

        Eigen::FullPivHouseholderQR<Eigen::MatrixXd> QR = J.fullPivHouseholderQr();
        delta = QR.solve(r) * std::pow(2.0, (double) armijo_int);
        long iter_rank = QR.rank();

        if (iter_rank < max_rank) {
            max_rank = iter_rank;
            std::cerr << "Rank = " << iter_rank << ", Equations = " << NumEquations << ", Variables = " << NumVariables << ", Iteration = " << j
                      << ". The jacobian rank is less than the number of equations. Some sketch elements may be over-constrained." << std::endl;

            // Randomly perturb the delta vector to prevent stagnation at inflection point
            std::random_device rdev;
            std::mt19937 rgen(rdev());
            auto rdist = std::uniform_real_distribution<>(-FLT_EPSILON * scale, FLT_EPSILON * scale);
            for (size_t i = 0; i != Variables.size(); ++i) {
                delta(Variables[i]->get_index()) += rdist(rgen);
            }
        }

        // Update Jacobian, residual, and Armijo step-size
        double iter_norm{DBL_MAX};
        int armijo_dir{0};
        while (iter_norm > residual_norm && residual_norm > std::sqrt(NumEquations + 1) * DBL_EPSILON) {
            for (size_t i = 0; i < Variables.size(); ++i) {
                Variables[i]->update(delta);
            }

            J.setZero();

            for (size_t i = 0; i != Verticies.size(); ++i) {
                Verticies[i]->update(J, r);
            }

            for (size_t i = 0; i != Curves.size(); ++i) {
                Curves[i]->update(J, r);
            }

            for (size_t i = 0; i != Constraints.size(); ++i) {
                Constraints[i]->update(J, r);
            }

            iter_norm = r.norm() / scale;

            if (armijo_dir == 0) {
                if (iter_norm < residual_norm) {
                    armijo_dir = +1;
                } else {
                    delta = -delta;
                    armijo_dir = -1;
                }
            }

            if (armijo_dir == 1) {
                if (iter_norm < residual_norm) {
                    if (armijo_int == 0) {
                        residual_norm = iter_norm;
                        break;
                    } else {
                        delta /= 2.0;
                        armijo_int += 1;
                    }
                } else {
                    delta *= -2.0;
                    armijo_dir = -1;
                }
            }

            if (armijo_dir == -1) {
                if (iter_norm > residual_norm) {
                    if (armijo_int == -(8 * sizeof(double))) {
                        std::cerr << "Residual norm was not reduced for Armijo step-sizes greater than 2^-(8*sizeof(double))" << std::endl;
                        residual_norm = iter_norm;
                        break;
                    } else {
                        delta /= 2.0;
                        armijo_int -= 1;
                    }
                } else {
                    residual_norm = iter_norm;
                    break;
                }
            }
        }
    }

    return residual_norm;
}

bool Sketch::build() {
    Constellation c = Constellation(this);

    Boundary = c.boundary();
    bool success = (Boundary->size() > 0);

    Contours = c.contours();
    success = success && (Contours.size() > 0);

    return success;
}

void Sketch::add_element(std::shared_ptr<Vertex> v) {
    add_element(Verticies, v);
}

void Sketch::add_element(std::shared_ptr<Curve> c) {
    add_element(Curves, c);
}

void Sketch::add_element(std::shared_ptr<Constraint> c) {
    add_element(Constraints, c);
}

void Sketch::add_element(std::shared_ptr<Pattern> p) {
    add_element(Patterns, p);
    p->register_elements(this);
}

template<>
void Sketch::save_as<SaveMethod::Rasterize>(std::string path, std::string file_name) const {
    /*
        This is a stub for visualization
    */

    if (!boost::filesystem::exists(path)) {
        boost::filesystem::create_directories(path);
    }

    std::fstream fs;
    fs.open(path + file_name + ".oesk", std::fstream::out);

    for (size_t i = 0; i < Curves.size(); ++i) {
        if (!Curves[i]->for_construction()) {
            for (size_t j = 0; j <= 10; ++j) {
                double2 v = Curves[i]->point(j / 10.0);
                fs << v.X << ',' << v.Y << '\n';
            }
            fs << "NaN" << ',' << "NaN" << '\n';
        }
    }

    fs.close();
}