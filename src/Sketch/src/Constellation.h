#ifndef CONSTELLATION_H
#define CONSTELLATION_H

#include "Sketch.h"

struct Branch {
	const Curve* Path;
	bool Orientation;
	double Angle;
};

class Star {
public:
	// Constructors
	Star() {};
	Star(const Vertex *v, const Sketch *s);

	const Vertex* vertex() const { return StarVertex; };

	const size_t size() const { return Branches.size(); };

	// Iterators
	std::list<Branch>::iterator begin() { return Branches.begin(); };
	std::list<Branch>::iterator end() { return Branches.end(); };

	std::list<Branch>::const_iterator begin() const { return Branches.begin(); };
	std::list<Branch>::const_iterator end() const { return Branches.end(); };

	std::list<Branch>::iterator next(std::list<Branch>::iterator iter) { return (iter == (--end()) ? begin() : ++iter); };
	std::list<Branch>::iterator prev(std::list<Branch>::iterator iter) { return (iter == begin() ? (--end()) : --iter); };

	std::list<Branch>::const_iterator next(std::list<Branch>::const_iterator iter) const { return (iter == (--end()) ? begin() : ++iter); };
	std::list<Branch>::const_iterator prev(std::list<Branch>::const_iterator iter) const { return (iter == begin() ? (--end()) : --iter); };

	// Curve pointer methods
	const Curve* next(const Curve* c) const;
	const Curve* prev(const Curve* c) const;

	void pop(const Curve* c);

private:
	const Vertex *StarVertex;
	std::list<Branch> Branches;
};

class Constellation {
public:
	// Constructors
	Constellation() {};
	Constellation(const Sketch *s);

	size_t size() { return Stars.size(); };
	
	bool contours(std::vector<Contour*> &contours);
	bool boundary(Contour* c);

private:
	std::list<Star> Stars;

	void pop(const Curve* c = nullptr);

	bool twin(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out);

	void supremum(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out);

	bool find_closed_contour(std::vector<const Curve*> &curves, std::vector<bool> &orientation);
};

#endif