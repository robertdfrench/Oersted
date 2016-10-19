#ifndef OERSTED_CONTOUR_H
#define OERSTED_CONTOUR_H

#include "Sketch.h"

class Contour {
public:
    // Constructors
    Contour() : Curves(std::vector<const Curve *>()), Orientation(std::vector<bool>()) {};

    Contour(const std::vector<const Curve *> &c);

    Contour(const std::vector<const Curve *> &c, const std::vector<bool> &dir);

    //Public Member Functions
    const Vertex *vertex(size_t i) const { return (Orientation[i] ? Curves[i]->start() : Curves[i]->end()); };

    const Curve *curve(size_t i) const { return Curves[i]; };

    const bool orientation(size_t i) const { return Orientation[i]; };

    const size_t size() const { return Curves.size(); };

    bool initialize(const std::vector<const Curve *> &c, const std::vector<bool> &dir);

    auto begin() { return Curves.begin(); };

    auto end() { return Curves.end(); };

    auto begin() const { return Curves.begin(); };

    auto end() const { return Curves.end(); };

    bool operator==(const Contour &rhs) const;

    double area() const;

private:
    std::vector<const Curve *> Curves;
    std::vector<bool> Orientation;
};

#endif //OERSTED_CONTOUR_H