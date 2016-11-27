#ifndef OERSTED_SKETCH_H
#define OERSTED_SKETCH_H

#include <cfloat>
#include <fstream>
#include <list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>

#include "Eigen"

class Variable;
class Vertex;
class Curve;
class Constraint;
class Pattern;
class Contour;

enum class SaveMethod {
    Rasterize
};

class Sketch {
public:
    Sketch() : NumVariables(0), NumEquations(0), Boundary() {};

    size_t size() const { return (Verticies.size() + Curves.size() + Constraints.size()); };

    size_t size_constraints() const { return Constraints.size(); };

    size_t size_contours() const { return Contours.size(); };

    size_t size_curves() const { return Curves.size(); };

    size_t size_verticies() const { return Verticies.size(); };

    bool build();

    double solve();

    void add_element(std::shared_ptr<Constraint> c);

    void add_element(std::shared_ptr<Curve> c);

    void add_element(std::shared_ptr<Pattern> p);

    void add_element(std::shared_ptr<Vertex> v);

    void add_parameter(std::shared_ptr<Variable> v) { add_parameter(Variables, v); };

    template<SaveMethod S>
    void save_as(std::string path, std::string file_name) const;

    template<class T, class...ArgT>
    std::shared_ptr<T> new_element(ArgT &&... args);

    std::shared_ptr<Curve const> curve(size_t i) const { return Curves[i]; };

    std::shared_ptr<Constraint const> constraint(size_t i) const { return Constraints[i]; };

    std::shared_ptr<Contour const> contour(size_t i) const { return Contours[i]; };

    std::shared_ptr<Contour const> boundary() const { return Boundary; };

    std::shared_ptr<Variable const> variable(size_t i) const { return Variables[i]; };

    std::shared_ptr<Vertex const> vertex(size_t i) const { return Verticies[i]; };

    std::vector<std::shared_ptr<Curve const>> const curves() { return std::vector<std::shared_ptr<Curve const>>{Curves.begin(),Curves.end()};};

private:
    std::vector<std::shared_ptr<Variable>> Variables;

    std::vector<std::shared_ptr<Vertex>> Verticies;

    std::vector<std::shared_ptr<Curve>> Curves;

    std::vector<std::shared_ptr<Constraint>> Constraints;

    std::vector<std::shared_ptr<Pattern>> Patterns;

    std::vector<std::shared_ptr<Contour>> Contours;

    std::shared_ptr<Contour> Boundary;

    size_t NumVariables;

    size_t NumEquations;

    template<class T>
    void add_element(std::vector<std::shared_ptr<T>> &elements, std::shared_ptr<T> e);

    template<class T>
    void add_parameter(std::vector<std::shared_ptr<T>> &parameters, std::shared_ptr<T> p);
};


template<class T>
void Sketch::add_element(std::vector<std::shared_ptr<T>> &elements, std::shared_ptr<T> e) {
    bool unique = true;
    for (size_t i = 0; i != elements.size() && unique; ++i) {
        unique = unique && (elements[i] != e);
    }

    if (unique) {
        NumEquations += e->set_equation_index(NumEquations);
        elements.push_back(e);
    }

    e->register_parameters(this);
}

template<class T>
void Sketch::add_parameter(std::vector<std::shared_ptr<T>> &parameters, std::shared_ptr<T> p) {
    bool unique = true;
    for (size_t i = 0; i != parameters.size() && unique; ++i) {
        unique = (parameters[i] != p);
    }

    if (unique) {
        NumVariables += p->set_index(NumVariables);
        parameters.push_back(p);
    }
}


template<class T, class...ArgT>
std::shared_ptr<T> Sketch::new_element(ArgT &&... args) {
    std::shared_ptr<T> e = std::make_shared<T>(std::forward<ArgT>(args)...);
    add_element(e);
    return e;
};

#endif //OERSTED_SKETCH_H