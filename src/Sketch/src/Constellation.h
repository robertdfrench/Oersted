#ifndef OERSTED_CONSTELLATION_H
#define OERSTED_CONSTELLATION_H

#include "Sketch.h"

class Constellation {
public:
    // Constructors
    Constellation() {};

    Constellation(const Sketch *s);

    size_t size() { return Stars.size(); };

    bool contours(std::vector<Contour *> &contours);

    bool boundary(Contour *c);

private:
    std::list<Star> Stars;

    void pop(const Curve *c = nullptr);

    bool twin(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out);

    void supremum(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out);

    bool find_closed_contour(std::vector<const Curve *> &curves, std::vector<bool> &orientation);
};

#endif //OERSTED_CONSTELLATION_H