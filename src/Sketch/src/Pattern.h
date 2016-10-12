#ifndef PATTERN_H
#define PATTERN_H

#include "Sketch.h"

class Pattern : public SketchElement {
public:
	void register_elements(Sketch *s);
	
	void register_parameters(Sketch *s) override {};
	size_t set_equation_index(size_t i) override { EquationIndex = i; return 0; };
	void update(Eigen::MatrixXd &J, Eigen::VectorXd &r) override {};

protected:
	std::vector<const Curve*> Input;

	std::vector<Vertex*> Verticies;
	std::vector<Curve*> Curves;
	std::vector<Constraint*> Constraints;
};

class MirrorCopy : public Pattern {
public:
	MirrorCopy(std::vector<const Curve*> &input, LineSegment *l);

private:
	LineSegment *SymmetryLine;
};

class RotateCopy : public Pattern {
public:
	RotateCopy(std::vector<const Curve*> &input, Vertex* center, double angle, size_t copies);

private:
	Vertex* Center;
	double Angle;
	size_t Copies;
};

#endif