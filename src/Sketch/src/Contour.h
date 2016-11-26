#ifndef OERSTED_CONTOUR_H
#define OERSTED_CONTOUR_H

#include <memory>
#include <vector>

#include "Curve.h"

class Contour {
public:
    // Constructors
    Contour() : Curves(std::vector<std::shared_ptr<Curve>>()), Orientation(std::vector<bool>()) {};

    Contour(std::vector<std::shared_ptr<Curve>> const &c);

    Contour(std::vector<std::shared_ptr<Curve>> const &c, std::vector<bool> const &dir);

    //Public Member Functions
    std::shared_ptr<Vertex> vertex(size_t i) const { return (Orientation[i] ? Curves[i]->start() : Curves[i]->end()); };

    std::shared_ptr<Curve> curve(size_t i) const { return Curves[i]; };

    const bool orientation(size_t i) const { return Orientation[i]; };

    const size_t size() const { return Curves.size(); };

    bool initialize(std::vector<std::shared_ptr<Curve>> const &c, std::vector<bool> const &dir);

    auto begin() { return Curves.begin(); };

    auto end() { return Curves.end(); };

    auto begin() const { return Curves.begin(); };

    auto end() const { return Curves.end(); };

    bool operator==(const Contour &rhs) const;

    double area() const;

private:
    std::vector<std::shared_ptr<Curve>> Curves;
    std::vector<bool> Orientation;
};

#endif //OERSTED_CONTOUR_H