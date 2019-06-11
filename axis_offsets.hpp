#include <vector>
#include "particle.hpp"
#include "crystal.hpp"

class axis_offsets {
    crystal * crystalp;
    std::vector<particle*> particles;
    std::vector<double> offsets;
    vec3 unit_direction;
public:
    axis_offsets(crystal * crystalp, vec3 direction) : crystalp(crystalp), unit_direction(direction.unit()) {
    }

    void add_particle(particle * p) {
        particles.push_back(p);
    }

    void measure() {
        offsets.clear();
        for (const particle * p: particles) {
            vec3 diff = crystalp->space.difference(p->cell->center, p->pos);
            double proj = diff * unit_direction;
            offsets.push_back(proj);
        }
    }

    double sum() {
        double total = 0;
        for (const auto & offset : offsets) {
            total += offset;
        }
        return total;
    }

    void write(std::ostream & out) const {
        bool first = true;
        for (const auto & offset : offsets) {
            if (!first) {
                out << " ";
            }
            out << offset;
            first = false;
        }
        out << std::endl;
    }
};
