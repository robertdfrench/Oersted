#ifndef OERSTED_STAR_H
#define OERSTED_STAR_H

#include <list>
#include <memory>

class Branch;
class Curve;
class Vertex;
class Sketch;

class Star {
public:
    Star() {};

    Star(std::shared_ptr<Vertex> v, Sketch const *s);

    size_t size() const { return Branches.size(); };

    void pop(std::shared_ptr<Curve> c);

    std::list<Branch>::const_iterator begin() const { return Branches.begin(); };

    std::list<Branch>::const_iterator end() const { return Branches.end(); };

    std::list<Branch>::const_iterator next(std::list<Branch>::const_iterator iter) const { return (iter == (--end()) ? begin() : ++iter); };

    std::list<Branch>::const_iterator prev(std::list<Branch>::const_iterator iter) const { return (iter == begin() ? (--end()) : --iter); };

    std::list<Branch>::iterator begin() { return Branches.begin(); };

    std::list<Branch>::iterator end() { return Branches.end(); };

    std::list<Branch>::iterator next(std::list<Branch>::iterator iter) { return (iter == (--end()) ? begin() : ++iter); };

    std::list<Branch>::iterator prev(std::list<Branch>::iterator iter) { return (iter == begin() ? (--end()) : --iter); };

    std::shared_ptr<Curve> next(std::shared_ptr<Curve> c) const;

    std::shared_ptr<Curve> prev(std::shared_ptr<Curve> c) const;

    std::shared_ptr<Vertex> const vertex() const { return StarVertex; };

private:
    std::shared_ptr<Vertex> StarVertex;
    std::list<Branch> Branches;
};

#endif //OERSTED_STAR_H
