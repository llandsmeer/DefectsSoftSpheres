#ifndef CRYSTAL_HPP
#define CRYSTAL_HPP

#include <iostream>
#include <fstream>

#include "vec3.hpp"
#include "particle.hpp"
#include "lattice_definition.hpp"
#include "lattice_cell.hpp"
#include "matrix3.hpp"
#include "periodic_space.hpp"

class crystal {
    crystal(periodic_space space) : space(space) { }
    void init(int n1, int n2, int n3, double a, lattice_definition & unitcell, double nncell_cutoff) {
        for (int i1 = 0; i1 < n1; i1++) {
            for (int i2 = 0; i2 < n2; i2++) {
                for (int i3 = 0; i3 < n3; i3++) {
                    for (size_t i4 = 0; i4 < unitcell.basis_vectors.size(); i4++) {
                        vec3 n(i1, i2, i3);
                        vec3 nb = n + unitcell.basis_vectors[i4];
                        lattice_cell * cell = new lattice_cell();
                        particle * p = new particle();
                        cell->n = n;
                        cell->nb = nb;
                        cell->owner = this;
                        p->pos = cell->center = space.project(nb);
                        p->cell = cell;
                        p->owner = this;
                        cell->basis = i4;
                        cell->particles.push_back(p);
                        cells.push_back(cell);
                        particles.push_back(p);
                    }
                }
            }
        }
        for (auto a : cells) {
            for (const auto b : cells) {
                double d = space.distance(a->center, b->center);
                if (a == b || d > nncell_cutoff) continue;
                a->nearest_neighbours.push_back(b);
            }
        }
    }
public:
    enum potential_type { STAR, HERTZ };
    potential_type potential_type = HERTZ;
    periodic_space space;
    std::vector<lattice_cell*> cells;
    std::vector<particle*> particles;

    double potential_sigma;
    double potential_epsilon; /* or f for star */

    double wigner_seitz_constraint = true;

    double potential(double dist) const {
        if (potential_type == HERTZ) {
            return dist >= potential_sigma ? 0 :
                potential_epsilon*pow(1.-dist/potential_sigma, 5./2);
        } else if (potential_type == STAR) {
            /* Phys. Rev. Let. V82 N26 p. 5290 eq 1
             * kBT is also in prefactor there but if we set
             * beta to 1 in the simulation its not needed */
            double f = potential_epsilon;
            double sig = potential_sigma;
            double prefactor = 5./18 * pow(f, 3./2);
            double _1_1psf2 = 1./(1.+sqrt(f)/2);
            if (dist <= sig) {
                return prefactor * (-std::log(dist/sig) + _1_1psf2);
            } else {
                return prefactor * (sig*_1_1psf2 * exp(-sqrt(f)*(dist-sig)/(2*sig))/dist);
            }
        }
        assert (false);
        return 0;
    }

    static crystal * bcc(int n1=4, int n2=-1, int n3=-1, double a=3) {
        if (n2 == -1) n2 = n1;
        if (n3 == -1) n3 = n2;
        auto unitcell = lattice_definition::body_centered_cubic(a);
        crystal * ret = new crystal(periodic_space(
            matrix3::from_cols(unitcell.p1(), unitcell.p2(), unitcell.p3()), vec3(n1, n2, n3)));
        ret->init(n1, n2, n3, a, unitcell, a*2); // 1.1
        return ret;
    }

    static crystal * hexagonal(int n1=4, int n2=-1, int n3=-1, double a=3) {
        if (n2 == -1) n2 = n1;
        if (n3 == -1) n3 = n2;
        auto unitcell = lattice_definition::hexagonal(a, 0.84);
        crystal * ret = new crystal(periodic_space(
            matrix3::from_cols(unitcell.p1(), unitcell.p2(), unitcell.p3()), vec3(n1, n2, n3)));
        ret->init(n1, n2, n3, a, unitcell, a*2);
        return ret;
    }

    double two_particle_energy(const particle * p1, const particle * p2,
            vec3 sh1=vec3(), vec3 sh2=vec3()) const {
        assert(p1 != p2);
        double energy_ws = 0;
#ifdef NDEBUG
        bool both = false;
#else
        bool both = wigner_seitz_constraint;
#endif
        if (both || wigner_seitz_constraint) {
            vec3 image1 = space.clip(p1->pos + sh1);
            for (const lattice_cell * nn : p1->cell->nearest_neighbours) {
                for (const particle * p : nn->particles) {
                    assert(p != p1 && nn != p1->cell);
                    if (p == p2) continue;
                    double dist = space.distance(image1, p->pos);
                    energy_ws += potential(dist);
                }
            }
            if (p1->cell->particles.size() > 1) {
                for (const particle * p : p1->cell->particles) {
                    if (p1 == p || p2 == p) continue;
                    double dist = space.distance(image1, p->pos);
                    energy_ws += potential(dist);
                }
            }
            vec3 image2 = space.clip(p2->pos + sh2);
            for (const lattice_cell * nn : p2->cell->nearest_neighbours) {
                for (const particle * p : nn->particles) {
                    assert(p != p2 && nn != p2->cell);
                    if (p == p1) continue;
                    double dist = space.distance(image2, p->pos);
                    energy_ws += potential(dist);
                }
            }
            if (p2->cell->particles.size() > 1) {
                for (const particle * p : p2->cell->particles) {
                    if (p2 == p || p1 == p) continue;
                    double dist = space.distance(image2, p->pos);
                    energy_ws += potential(dist);
                }
            }
            energy_ws += 2*potential(space.distance(image1, image2));
#ifdef NDEBUG
            return energy_ws;
#endif
        }
        double energy_all = 0;
        if (both || !wigner_seitz_constraint) {
            vec3 image1 = space.clip(p1->pos + sh1);
            vec3 image2 = space.clip(p2->pos + sh2);
            for (const particle * p : particles) {
                if (p == p1 || p == p2) continue;
                double dist1 = space.distance(image1, p->pos);
                energy_all += potential(dist1);
                double dist2 = space.distance(image2, p->pos);
                energy_all += potential(dist2);
            }
            energy_all += 2*potential(space.distance(image1, image2));
#ifdef NDEBUG
            return energy_all;
#endif
        }
        if (!(abs(energy_ws - energy_all) < 1e-3)) {
            std::cout << "two_particle_energy() mismatch\n";
            PRINT_VAR(p1->cell->nb);
            PRINT_VAR(p2->cell->nb);
            PRINT_VAR(energy_ws);
            PRINT_VAR(energy_all);
        }
        assert(abs(energy_ws - energy_all) < 1e-3);
        return energy_all;
    }

    /* boring functions */

    lattice_cell * get_cell(int n1, int n2, int n3, int n4) {
        vec3 n(n1, n2, n3);
        for (lattice_cell * cell : cells) {
            if (cell->n == n && cell->basis == n4) {
                return cell; }
        }
        throw "cell not found";
    }

    double density() const {
        return particles.size() / space.volume();
    }

    void write(const std::string & filename) const {
        std::ofstream stream(filename);
        write(stream);
        stream.close();
    }

    void write(std::ostream & out) const {
        out << particles.size() << "\n";
        auto p1 = space.p1();
        auto p2 = space.p2();
        auto p3 = space.p3();
        out << "0 0 0" << "\n";
        out << p1.x << " " << p1.y << " " << p1.z << "\n";
        out << p2.x << " " << p2.y << " " << p2.z << "\n";
        out << p3.x << " " << p3.y << " " << p3.z << "\n";
        for (const particle * p: particles) {
            out << p->pos.x << " "
                << p->pos.y << " "
                << p->pos.z << " "
                << p->size << " "
                << p->color << " "
                << "\n"; 
        }
        out << std::flush;
    }

    void log(int iter, std::ostream & out) const {
        for (const particle * p: particles) {
            out << iter << " "
                << p->pos.x << " "
                << p->pos.y << " "
                << p->pos.z << " "
                << p->cell->center.x << " "
                << p->cell->center.y << " "
                << p->cell->center.z << " "
                << p->cell->n.x << " "
                << p->cell->n.y << " "
                << p->cell->n.z << " "
                << p->cell->basis << " "
                << "\n"; 
        }
        out << std::flush;
    }

};

// stupid c++

particle * lattice_cell::interstitial(vec3 offset) {
    assert(contains(center + offset));
    for (const particle * p : particles) {
        contains(p->pos - offset / particles.size());
    }
    for (particle * p : particles) {
        p->pos = owner->space.clip(p->pos - offset / particles.size());
    }
    particle * p = new particle();
    p->pos = owner->space.clip(center + offset);
    p->cell = this;
    p->owner = owner;
    owner->particles.push_back(p);
    particles.push_back(p);
    return p;
}

void lattice_cell::vacancy(int basis) {
    assert((int)particles.size() >= basis);
    particle * p = particles[basis];
    for (auto it = owner->particles.begin(); it != owner->particles.end();) {
        if ((*it) == p) {
            it = owner->particles.erase(it);
        } else {
            ++it;
        }
    }
    particles.erase(particles.begin() + basis);
}

double particle::energy(vec3 shift) {
#ifdef NDEBUG
        bool both = false;
#else
        bool both = owner->wigner_seitz_constraint;
#endif
    assert(cell->contains(pos));
    double energy_ws = 0;
    if (both || owner->wigner_seitz_constraint) {
        for (const lattice_cell * nn : cell->nearest_neighbours) {
            vec3 image = owner->space.clip(pos + shift);
            for (const particle * p : nn->particles) {
                assert(p != this);
                double dist = owner->space.distance(image, p->pos);
                energy_ws += owner->potential(dist);
            }
        }
        if (cell->particles.size() > 1) {
            vec3 image = owner->space.clip(pos + shift);
            for (const particle * p : cell->particles) {
                if (p == this) continue;
                double dist = owner->space.distance(image, p->pos);
                energy_ws += owner->potential(dist);
            }
        }
#ifdef NDEBUG
        return energy_ws;
#endif
    }
    double energy_all = 0;
    if (both || !owner->wigner_seitz_constraint) {
        vec3 image = owner->space.clip(pos + shift);
        for (const particle * p : owner->particles) {
            if (p == this) continue;
            double dist = owner->space.distance(image, p->pos);
            energy_all += owner->potential(dist);
        }
#ifdef NDEBUG
        return energy_all;
#endif
    }
    assert(energy_ws == energy_all);
    return energy_all;
}

bool lattice_cell::contains(const vec3 & pos) const {
    /* wigner seitz constraint */
    if (!owner->wigner_seitz_constraint) return true;
    auto d1 = owner->space.distance(center, pos);
    for (const lattice_cell * nn : nearest_neighbours) {
        auto d2 = owner->space.distance(nn->center, pos);
        if (d2 < d1) {
            return false;
        }
    }
    return true;
}

#endif
