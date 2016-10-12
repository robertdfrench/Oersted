#include "Sketch.hpp"

Star::Star(const Vertex *v, const Sketch *s) {
	StarVertex = v;

	// Extract Curves
	std::vector<const Curve*> curves;
	std::vector<bool> orientation;
	for (size_t i = 0;i < s->size_curves();++i) {
		const Curve* c = s->curve(i);
		if (!c->ForConstruction) {
			if (StarVertex == c->start()) {
				curves.push_back(c);
				orientation.push_back(true);
			}
			else if (StarVertex == c->end()) {
				curves.push_back(c);
				orientation.push_back(false);
			}
		}
	}

	// Calculate angle tangent to star point
	std::vector<double> angle;
	std::vector<double> d_angle;
	for (size_t i = 0;i < curves.size();++i) {
		if (orientation[i]) {
			angle.push_back(curves[i]->a(0.0, true));
			d_angle.push_back(curves[i]->da(0.0, true));
		}
		else {
			angle.push_back(curves[i]->a(1.0, false));
			d_angle.push_back(curves[i]->da(1.0, false));
		}
	}

	// Sort by tangent angle
	std::vector<size_t> index;
	index.resize(angle.size());
	std::iota(index.begin(), index.end(), 0);
	
	double tol = M_PI * sqrt(DBL_EPSILON);
	std::sort(index.begin(), index.end(),
		[&](size_t i, size_t j) {	
			if (abs(angle[i] - angle[j]) < tol || (abs(abs(angle[i]) - M_PI) < tol && abs(abs(angle[j]) - M_PI) < tol)) {
				return (d_angle[i] > d_angle[j]);
			}
			else {
				return (angle[i] > angle[j]);
			}
		}
	);

	// Store end point angle difference
	double x = StarVertex->x();
	double y = StarVertex->y();
	for (size_t i = 0;i != curves.size();++i) {
		const Vertex *v;
		if (orientation[i]) {
			v = curves[i]->end();
		}
		else {
			v = curves[i]->start();
		}

		double dx = v->x() - x;
		double dy = v->y() - y;

		angle[i] = atan2(dy, dx);
	}

	for (size_t i = 0;i != index.size();++i) {
		size_t j = index[i];
		size_t k = index[(i + 1) % index.size()];

		Branches.push_back(Branch{ curves[j], orientation[j], angle[j] - angle[k] });

		// Handle branch cut
		if (Branches.back().Angle < 0.0) {
			Branches.back().Angle += 2.0 * M_PI;
		}
	}

	// Handle self-closed Star with two branches having shared endpoints
	if (index.size() == 2 && Branches.back().Angle == 0.0) {
		Branches.back().Angle += 2.0 * M_PI;
	}
	
}

const Curve* Star::next(const Curve* c) const {
	for (auto b = begin(); b != end(); ++b) {
		if (b->Path == c) {
			return next(b)->Path;
		}
	}

	return nullptr;
}

const Curve* Star::prev(const Curve* c) const {
	for (auto b = begin(); b != end(); ++b) {
		if (b->Path == c) {
			return prev(b)->Path;
		}
	}

	return nullptr;
}

void Star::pop(const Curve* c) {
	for (auto b = begin();b != end();++b) {
		if (b->Path == c) {
			prev(b)->Angle += b->Angle;

			Branches.erase(b);

			break;
		}
	}
}

Constellation::Constellation(const Sketch *s) {
	for (size_t i = 0;i != s->size_verticies();++i) {
		Stars.push_back(Star{ s->vertex(i),s });
	}
	pop();
}

bool Constellation::twin(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out) {
	const Curve* c = b_out->Path;

	for (auto s = Stars.begin();s != Stars.end();++s) {
		if (s != s_out) {
			for (auto b = s->begin();b != s->end();++b) {
				if (b->Path == b_out->Path) {
					b_out = b;
					s_out = s;

					return true;
				}
			}
		}
	}

	return false;
}

void Constellation::supremum(std::list<Star>::iterator &s_out, std::list<Branch>::iterator &b_out) {
	double sup = 0.0;
	double ang;

	for (auto s = Stars.begin();s != Stars.end();++s) {
		for (auto b = s->begin();b != s->end();++b) {
			double sup_ij = b->Path->supremum();
			if (sup < sup_ij || (sup == sup_ij && b->Angle < ang)) {
				sup = sup_ij;
				ang = b->Angle;

				s_out = s;
				b_out = b;
			}
		}
	}
}

void Constellation::pop(const Curve* c) {
	bool iterate = true;

	while (iterate) {
		iterate = false;

		// remove curve c
		for (auto s = Stars.begin();s != Stars.end();++s) {
			s->pop(c);
		}

		auto s = Stars.begin();
		while (s!= Stars.end()) {
			if (s->size() == 0) {
				s = Stars.erase(s);
			}
			else {
				++s;
			}
		}

		// pop curves from stars which only have one branch
		for (auto s = Stars.begin();s != Stars.end();++s) {
			if (s->size() == 1) {
				c = s->begin()->Path;
				iterate = true;
				break;
			}
		}
	}
}

bool Constellation::boundary(Contour* c) {
	std::vector<const Curve*> curves;
	std::vector<bool> orientation;

	// Base of induction
	auto s{ Stars.begin() };
	auto b{ s->begin() };

	supremum(s, b);

	double angle{ 0.0 };
	b = s->prev(b);

	curves.push_back(b->Path);
	orientation.push_back(b->Orientation);

	while (true) {
		bool has_twin = twin(s, b);

		if (!has_twin) {
			return false;
		}
		else {
			b = s->prev(b);
			angle += 2.0 * M_PI - b->Angle;
		}

		if (b->Path == curves.back()) { //excludes self-closed contours
			return false;
		}
		else if (b->Path == curves[0]) {
			double tol = (FLT_EPSILON * M_PI * (curves.size() - 1));
			double expected = (curves.size() - 2) * M_PI;
			double diff = abs(angle - expected);

			if (diff <= tol) {
				c->initialize(curves, orientation);
				return true;
			}
			else {
				return false;
			}
		}
		else {
			curves.push_back(b->Path);
			orientation.push_back(b->Orientation);
		}
	}
}

bool Constellation::contours(std::vector<Contour*> &contours) {
	std::vector<const Curve*> contour_curves;
	std::vector<bool> orientation;

	while (size() > 0) {
		bool success = find_closed_contour(contour_curves, orientation);
		if (success) {
			contours.push_back(new Contour(contour_curves, orientation));
		}
		else {
			return false;
		}
	}

	return true;
}

bool Constellation::find_closed_contour(std::vector<const Curve*> &curves, std::vector<bool> &orientation) {
	curves.resize(0);
	orientation.resize(0);

	// Base of induction
	auto s{ Stars.begin() };
	auto b{ s->begin() };

	supremum(s, b);

	double angle{0.0};
	b = s->next(b);

	curves.push_back(b->Path);
	orientation.push_back(b->Orientation);

	while (true) {
		bool has_twin = twin(s, b);

		if (!has_twin) {
			return false;
		}
		else {
			angle += b->Angle;
			b = s->next(b);
		}

		if (b->Path == curves.back()) { //excludes self-closed contours
			return false;
		}
		else if (b->Path == curves[0]) {
			double tol = (FLT_EPSILON * M_PI * (curves.size() - 1));
			double expected = (curves.size() - 2) * M_PI;
			double diff = abs(angle - expected);

			if (diff <= tol) {
				pop(curves.back());
				return true;
			}
			else {
				return false;
			}
		}
		else {
			curves.push_back(b->Path);
			orientation.push_back(b->Orientation);
		}
	}
}