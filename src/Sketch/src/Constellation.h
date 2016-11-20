#ifndef OERSTED_CONSTELLATION_H
#define OERSTED_CONSTELLATION_H

#include "Sketch.h"

class Constellation {
public:
    // Constructors
    Constellation() {};

    Constellation(Sketch const *s);

    size_t size() { return Stars.size(); };

    std::vector<std::shared_ptr<Contour>> contours();

    std::shared_ptr<Contour> boundary();

private:
    std::list<Star> Stars;

    void pop(std::shared_ptr<Curve> c = std::shared_ptr<Curve>());

    bool twin(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out);

    void supremum(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out);

    bool find_closed_contour(std::vector<std::shared_ptr<Curve>> &curves, std::vector<bool> &orientation);
};

#endif //OERSTED_CONSTELLATION_H