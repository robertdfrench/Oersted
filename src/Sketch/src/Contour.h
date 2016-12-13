#ifndef OERSTED_CONTOUR_H
#define OERSTED_CONTOUR_H

#include <memory>
#include <vector>

#include "Curve.h"

class Contour {
public:
    Contour() : Curves(std::vector<std::shared_ptr<Curve const>>()), Orientation(std::vector<bool>()) {};

    Contour(std::vector<std::shared_ptr<Curve const>> const &c);

    Contour(std::vector<std::shared_ptr<Curve const>> const &c, std::vector<bool> const &dir);

    size_t size() const { return Curves.size(); };

    auto begin() { return Curves.begin(); };

    auto end() { return Curves.end(); };

    auto begin() const { return Curves.begin(); };

    auto end() const { return Curves.end(); };

    bool initialize(std::vector<std::shared_ptr<Curve const>> const &c, std::vector<bool> const &dir);

    bool orientation(size_t i) const { return Orientation[i]; };

    bool operator==(Contour const &rhs) const;

    double area() const;

    std::shared_ptr<Vertex const> vertex(size_t i) const { return (Orientation[i] ? Curves[i]->start() : Curves[i]->end()); };

    std::shared_ptr<Curve const> curve(size_t i) const { return Curves[i]; };

private:
    std::vector<std::shared_ptr<Curve const>> Curves;

    std::vector<bool> Orientation;
};

#endif //OERSTED_CONTOUR_H