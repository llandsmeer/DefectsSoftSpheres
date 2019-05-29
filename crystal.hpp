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
public:
    periodic_space space;
    std::vector<lattice_cell*> cells;
    std::vector<particle*> particles;

    double potential_sigma;
    double potential_epsilon;

    double wigner_seitz_constraint = true;

    double potential(double dist) const {
        return dist >= potential_sigma ? 0 :
            potential_epsilon*pow(1.-dist/potential_sigma, 5./2);
    }

    static crystal * bcc(int n1=4, int n2=-1, int n3=-1, double a=3) {
        if (n2 == -1) n2 = n1;
        if (n3 == -1) n3 = n2;
        auto unitcell = lattice_definition::body_centered_cubic(a);
        crystal * ret = new crystal(periodic_space(
            matrix3::from_cols(
                unitcell.p1(),
                unitcell.p2(),
                unitcell.p3()),
            vec3(n1, n2, n3)));
        for (int i1 = 0; i1 < n1; i1++) {
            for (int i2 = 0; i2 < n2; i2++) {
                for (int i3 = 0; i3 < n3; i3++) {
                    for (int i4 = 0; i4 < unitcell.basis_vectors.size(); i4++) {
                        vec3 n(i1, i2, i3);
                        vec3 nb = n + unitcell.basis_vectors[i4];
                        lattice_cell * cell = new lattice_cell();
                        particle * p = new particle();
                        cell->n = n;
                        cell->nb = nb;
                        cell->owner = ret;
                        p->pos = cell->center = ret->space.project(nb);
                        p->cell = cell;
                        p->owner = ret;
                        cell->basis = i4;
                        cell->particles.push_back(p);
                        ret->cells.push_back(cell);
                        ret->particles.push_back(p);
                    }
                }
            }
        }
        auto nncell_cutoff = a*1.1;
        for (auto a : ret->cells) {
            for (const auto b : ret->cells) {
                double d = ret->space.distance(a->center, b->center);
                if (a == b || d > nncell_cutoff) continue;
                a->nearest_neighbours.push_back(b);
            }
        }
        return ret;
    }

    double two_particle_energy(const particle * p1, const particle * p2,
            vec3 sh1=vec3(), vec3 sh2=vec3()) const {
        assert(p1 != p2);
        double energy = 0;
        if (wigner_seitz_constraint) {
            vec3 image = space.clip(p1->pos + sh1);
            for (const lattice_cell * nn : p1->cell->nearest_neighbours) {
                for (const particle * p : nn->particles) {
                    assert(p != p1);
                    double dist = space.distance(image, p->pos);
                    energy += potential(dist);
                }
            }
            if (p1->cell->particles.size() > 1) {
                for (const particle * p : p1->cell->particles) {
                    if (p1 == p) continue;
                    double dist = space.distance(image, p1->pos);
                    energy += potential(dist);
                }
            }
            image = space.clip(p2->pos + sh2);
            for (const lattice_cell * nn : p2->cell->nearest_neighbours) {
                for (const particle * p : nn->particles) {
                    assert(p != p2);
                    double dist = space.distance(image, p->pos);
                    energy += potential(dist);
                }
            }
            if (p2->cell->particles.size() > 1) {
                for (const particle * p : p2->cell->particles) {
                    if (p2 == p) continue;
                    double dist = space.distance(image, p2->pos);
                    energy += potential(dist);
                }
            }
        } else {
            vec3 image1 = space.clip(p1->pos + sh1);
            vec3 image2 = space.clip(p2->pos + sh2);
            for (const particle * p : particles) {
                if (p != p1) {
                    double dist1 = space.distance(image1, p->pos);
                    energy += potential(dist1);
                }
                if (p != p2) {
                    double dist2 = space.distance(image2, p->pos);
                    energy += potential(dist2);
                }
            }
        }
        return energy;
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

};

// stupid c++

particle * lattice_cell::interstitial(vec3 offset) {
    assert(contains(center + offset));
    for (const particle * p : particles) {
        contains(p->pos + offset / particles.size());
    }
    particle * p = new particle();
    p->pos = center + offset;
    p->cell = this;
    p->owner = owner;
    owner->particles.push_back(p);
    particles.push_back(p);
    return p;
}

double particle::energy(vec3 shift) {
    assert(cell->contains(pos));
    double energy = 0;
    if (owner->wigner_seitz_constraint) {
        for (const lattice_cell * nn : cell->nearest_neighbours) {
            vec3 image = owner->space.clip(pos + shift);
            for (const particle * p : nn->particles) {
                assert(p != this);
                double dist = owner->space.distance(image, p->pos);
                energy += owner->potential(dist);
            }
        }
        if (cell->particles.size() > 1) {
            vec3 image = owner->space.clip(pos + shift);
            for (const particle * p : cell->particles) {
                if (p == this) continue;
                double dist = owner->space.distance(image, p->pos);
                energy += owner->potential(dist);
            }
        }
    } else {
        vec3 image = owner->space.clip(pos + shift);
        for (const particle * p : owner->particles) {
            if (p == this) continue;
            double dist = owner->space.distance(image, p->pos);
            energy += owner->potential(dist);
        }
    }
    return energy;
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
