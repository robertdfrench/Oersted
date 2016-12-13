#include "Contour.h"

Contour::Contour(std::vector<std::shared_ptr<Curve const>> const &c) {
    // TODO: Check for intersecting curves
    // TODO: Check for ccw orientation of entire contour

    Curves = std::vector<std::shared_ptr<Curve const>>();
    Curves.reserve(c.size());

    std::vector<std::shared_ptr<Vertex const>> start;
    start.reserve(c.size());

    std::vector<std::shared_ptr<Vertex const>> end;
    end.reserve(c.size());

    for (size_t i = 0; i < c.size(); i++) {
        Curves.push_back(c[i]);
        start.push_back(c[i]->start());
        end.push_back(c[i]->end());
    }

    Orientation = std::vector<bool>(c.size(), true);

    bool simple = true;
    double angle = 0.0;

    for (size_t i = 0; i < Curves.size() - 1; i++) {
        for (size_t j = i + 1; j < Curves.size(); j++) {
            if (end[i] == start[j]) {
                std::swap(Curves[i + 1], Curves[j]);
                std::swap(start[i + 1], start[j]);
                std::swap(end[i + 1], end[j]);
            } else if (end[i] == end[j]) {
                std::swap(Curves[i + 1], Curves[j]);
                std::swap(start[i + 1], end[j]);
                std::swap(end[i + 1], start[j]);
                std::swap(end[j], start[j]);

                std::swap(Orientation[i + 1], Orientation[j]);
                Orientation[i + 1] = !Orientation[i + 1];
            }
        }

        if (end[i] != start[i + 1]) {
            simple = false;
            break;
        }
    }

    bool closed = (start.front() == end.back());

    if (!simple) {
        throw std::exception(); // Curves do not form a simple contour
    } else if (!closed) {
        throw std::exception(); // Curves do not form a closed contour
    }
}

Contour::Contour(std::vector<std::shared_ptr<Curve const>> const &c, std::vector<bool> const &dir) {
    bool success = initialize(c, dir);

    if (!success) {
        throw std::exception(); // Curve are not correctly oriented or do not form a closed contour
    }
}

bool Contour::initialize(std::vector<std::shared_ptr<Curve const>> const &c, std::vector<bool> const &dir) {
    Curves = std::vector<std::shared_ptr<Curve const>>();
    Curves.reserve(c.size());

    Orientation = std::vector<bool>();
    Orientation.reserve(c.size());

    for (size_t i = 0; i != c.size(); ++i) {
        Curves.push_back(c[i]);
        Orientation.push_back(dir[i]);

        size_t j = (i + 1) % c.size();
        std::shared_ptr<Vertex const> vi = (dir[i] ? c[i]->end() : c[i]->start());
        std::shared_ptr<Vertex const> vj = (dir[j] ? c[j]->start() : c[j]->end());
        if (vi != vj) {
            return false;
        }
    }

    return true;
}

bool Contour::operator==(const Contour &rhs) const {
    if (size() != rhs.size()) {
        return false;
    }

    bool match = false;
    bool direction;
    size_t d;
    for (size_t i = 0; i < Curves.size(); ++i) {
        if (Curves[i] == rhs.Curves[0]) {
            match = true;
            direction = (Orientation[i] == rhs.Orientation[0]);
            d = i;
            break;
        }
    }

    if (match) {
        for (size_t j = 0; j < Curves.size(); ++j) {
            size_t i = (direction ? (d + j) : (d - j)) % Curves.size();

            if (Curves[i] != rhs.Curves[j]) {
                return false;
            }
        }

        return true;
    } else {
        return false;
    }
}

double Contour::area() const {
    double area{0.0};

    for (size_t i = 0; i != Curves.size(); ++i) {
        std::shared_ptr<Vertex const> v0 = Curves[i]->start();
        std::shared_ptr<Vertex const> v1 = Curves[i]->end();
        double da = Curves[i]->area() + v0->x() * v1->y() - v0->y() * v1->x();

        if (Orientation[i]) {
            area += da;
        } else {
            area -= da;
        }
    }

    return area;
}