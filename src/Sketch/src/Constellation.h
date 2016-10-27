#ifndef OERSTED_CONSTELLATION_H
#define OERSTED_CONSTELLATION_H

#include "Sketch.h"

class Constellation {
public:
    // Constructors
    Constellation() {};

    Constellation(const Sketch *s);

    size_t size() { return Stars.size(); };

    bool contours(vector<Contour *> &contours);

    bool boundary(Contour *c);

private:
    list<Star> Stars;

    void pop(const Curve *c = nullptr);

    bool twin(list<Star>::iterator &s_out, list<Branch>::iterator &b_out);

    void supremum(list<Star>::iterator &s_out, list<Branch>::iterator &b_out);

    bool find_closed_contour(vector<const Curve *> &curves, vector<bool> &orientation);
};

#endif //OERSTED_CONSTELLATION_H