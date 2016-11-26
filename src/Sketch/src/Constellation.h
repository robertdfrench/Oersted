#ifndef OERSTED_CONSTELLATION_H
#define OERSTED_CONSTELLATION_H

#include <list>
#include <vector>
#include <memory>

class Branch;
class Contour;
class Curve;
class Star;
class Sketch;

class Constellation {
public:
    Constellation() {};

    Constellation(Sketch const *s);

    size_t size() { return Stars.size(); };

    std::vector<std::shared_ptr<Contour>> contours();

    std::shared_ptr<Contour> boundary();

private:
    std::list<Star> Stars;

    bool find_closed_contour(std::vector<std::shared_ptr<Curve>> &curves, std::vector<bool> &orientation);

    bool twin(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out);

    void pop(std::shared_ptr<Curve> c = std::shared_ptr<Curve>());

    void supremum(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out);
};

#endif //OERSTED_CONSTELLATION_H