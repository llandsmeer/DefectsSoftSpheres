#ifndef AXIS_OFFSETS_HPP
#define AXIS_OFFSETS_HPP

#include <array>
#include <functional>
#include <vector>
#include <algorithm>
#include "particle.hpp"
#include "crystal.hpp"

class axis_offsets {
    crystal * crystalp;
    vec3 unit_direction;
public:
    int label;
    std::vector<double> offsets;
    struct pref {
        particle * p;
        int state; /* -1, 0 or +1, relative crowdion position after trace */
        pref(particle * p, int state) : p(p), state(state) { }
    };
    std::vector<pref> particles;
    axis_offsets(crystal * crystalp, vec3 direction) : crystalp(crystalp), unit_direction(direction.unit()) {
    }

    void add_particle(particle * p, int state=0) {
        particles.emplace_back(p, state);
    }

    void measure() {
        offsets.clear();
        for (auto & ref: particles) {
            vec3 diff = crystalp->space.difference(ref.p->cell->center, ref.p->pos);
            double proj = diff * unit_direction;
            offsets.push_back(proj);
        }
    }

    double sum() {
        double total = 0;
        for (const auto & offset : offsets) {
            total += abs(offset);
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

    void add_trace(lattice_cell * cell, double radius=0.5) {
        vec3 center = cell->center;
        assert(particles.size() == 0);
        vec3 probe = center;
        int niter = 0;
        for (;;) {
            if (niter > 10 && crystalp->space.distance(probe, center) < radius) {
                break;
            }
            for (particle * p : crystalp->particles) {
                if (crystalp->space.distance(probe, p->cell->center) < radius) {
                    bool already_added = false;
                    for (const auto & qref: particles) {
                        if (p == qref.p) {
                            already_added = true;
                            break;
                        }
                    }
                    if (!already_added) {
                        double d = crystalp->space.difference(center, p->cell->center) * unit_direction;
                        int state = d > 0 ? 1 : -1;
                        if (p->cell == cell) state = 0;
                        particles.emplace_back(p, state);
                    }
                }
            }
            probe = crystalp->space.clip(probe + radius * unit_direction);
            niter += 1;
        }
        std::sort(std::begin(particles), std::end(particles), [&](pref & a, pref & b) -> bool {
                double ad = crystalp->space.difference(center, a.p->cell->center) * unit_direction;
                double bd = crystalp->space.difference(center, b.p->cell->center) * unit_direction;
                return ad < bd;
        });
    }
};

#endif
