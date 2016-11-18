#ifndef OERSTED_STAR_H
#define OERSTED_STAR_H

#include "Sketch.h"

class Star {
public:
    // Constructors
    Star() {};

    Star(std::shared_ptr<Vertex> v, Sketch const *s);

    std::shared_ptr<Vertex> const vertex() const { return StarVertex; };

    const size_t size() const { return Branches.size(); };

    // Iterators
    std::list<Branch>::iterator begin() { return Branches.begin(); };

    std::list<Branch>::iterator end() { return Branches.end(); };

    std::list<Branch>::const_iterator begin() const { return Branches.begin(); };

    std::list<Branch>::const_iterator end() const { return Branches.end(); };

    std::list<Branch>::iterator next(std::list<Branch>::iterator iter) {
        return (iter == (--end()) ? begin() : ++iter);
    };

    std::list<Branch>::iterator prev(std::list<Branch>::iterator iter) {
        return (iter == begin() ? (--end()) : --iter);
    };

    std::list<Branch>::const_iterator next(std::list<Branch>::const_iterator iter) const {
        return (iter == (--end()) ? begin() : ++iter);
    };

    std::list<Branch>::const_iterator prev(std::list<Branch>::const_iterator iter) const {
        return (iter == begin() ? (--end()) : --iter);
    };

    // Curve pointer methods
    const Curve *next(const Curve *c) const;

    const Curve *prev(const Curve *c) const;

    void pop(const Curve *c);

private:
    std::shared_ptr<Vertex> StarVertex;
    std::list<Branch> Branches;
};

#endif //OERSTED_STAR_H
