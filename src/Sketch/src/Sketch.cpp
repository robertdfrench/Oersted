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
    // TODO: Use sparse matrix representation for Jacobian, (this is lower priority since typical uses will not produce very large matrices)
    // TODO: Try extracting sub matrix and solving (Certain constraints are known to eliminate variables (Radius, Fixation), if the matrix can be permuted so that a submatrix on the block diagonal is lower-triangular, then a subsystem of equations can be solved independently of the remaining equations. Moreover, if the jacobian can be permuted into a block lower triangular form, subsets of equations can be solved blockwise)
    // TODO: Use variable_feature_size() and equation_feature_size()

    // Matrix and vectors
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(NumEquations, NumVariables);
    Eigen::VectorXd r = Eigen::VectorXd::Zero(NumEquations);
    Eigen::VectorXd lfs = Eigen::VectorXd::Zero(NumEquations);

    // Tolerances and iteration limits
    double dbl_tol{std::sqrt(NumEquations + 1.0) * DBL_EPSILON};
    double flt_tol{std::sqrt(NumEquations + 1.0) * FLT_EPSILON};

    size_t double_bits = (CHAR_BIT * sizeof(double));

    // Initial residual calculation
    update_linearization(J, r);

    // Setup iteration
    double res_norm{r.norm()};
    double res_tol{characteristic_length() * dbl_tol};
    long prev_rank(std::min(NumEquations, NumVariables));
    for (size_t j = 0; res_norm > res_tol && j != double_bits; ++j) {
        Eigen::FullPivHouseholderQR<Eigen::MatrixXd> QR = J.fullPivHouseholderQr();
        Eigen::VectorXd delta = QR.solve(r);

        if (((J * delta) - r).norm() > res_tol + r.norm() * flt_tol) {
            std::cerr << "Linearized sketch equation could not be solved with an average error per equation less than FLT_EPSILON. Some elements may be over-constrained." << std::endl;
        }

        long iter_rank = QR.rank();
        if (iter_rank < prev_rank) {
            // perturb(delta, characteristic_length()); // perturb update in case of stagnation at inflection point
        }
        prev_rank = iter_rank;

        update_variables(delta);
        update_linearization(J, r);

        double iter_norm{r.norm()};
        double iter_tol{res_tol};
        if (iter_norm > res_norm) { // Line search
            delta *= -0.5;

            size_t line_search_iter{0};
            while (iter_norm > res_norm && iter_norm > iter_tol && line_search_iter != double_bits) {
                update_variables(delta);
                update_linearization(J, r);

                iter_norm = r.norm();
                iter_tol = dbl_tol * characteristic_length();
                if (iter_norm > res_norm) {
                    delta *= 0.5;
                    line_search_iter += 1;

                    if (line_search_iter == double_bits) {
                        std::cerr << "Residual norm was not reduced for any step-size greater than 2^-(CHAR_BIT * sizeof(double))" << std::endl;
                    }
                }
            }
        }
        res_norm = iter_norm;
        res_tol = iter_tol;
    }

    return res_norm;
}

void Sketch::update_linearization(Eigen::MatrixXd &J, Eigen::VectorXd &r) const {
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
}

void Sketch::update_variables(Eigen::VectorXd &delta) const {
    for (size_t i = 0; i != Variables.size(); ++i) {
        Variables[i]->update(delta);
    }
}

void Sketch::perturb(Eigen::VectorXd &delta, double scale) const {
    // Randomly perturb the delta vector to prevent stagnation at inflection point
    std::random_device rdev;
    std::mt19937 rgen(rdev());
    std::uniform_real_distribution<double> rdist = std::uniform_real_distribution<double>(-DBL_EPSILON * scale, DBL_EPSILON * scale);
    for (size_t i = 0; i != Variables.size(); ++i) {
        size_t j = Variables[i]->get_index();
        delta(j) += rdist(rgen) * scale;
    }
}

double Sketch::characteristic_length() const {
    double scale{0.0};

    if (Curves.size() > 0) {
        for (size_t i = 0; i != Curves.size(); ++i) {
            scale = std::max(scale, Curves[i]->length());
        }
    } else if (Verticies.size() > 0) {
        double xmax = DBL_MIN;
        double xmin = DBL_MAX;
        double ymax = DBL_MIN;
        double ymin = DBL_MAX;
        for (size_t i = 0; i != Verticies.size(); ++i) {
            xmax = std::max(xmax, Verticies[i]->x());
            xmin = std::min(xmin, Verticies[i]->x());

            ymax = std::max(ymax, Verticies[i]->y());
            ymin = std::max(ymin, Verticies[i]->y());
        }

        scale = std::max(xmax - xmin, ymax-ymin);
        if(scale == 0) {
            scale = std::max(xmax,ymax);
        }
    }

    return scale;
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