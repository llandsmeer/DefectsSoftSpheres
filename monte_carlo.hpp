#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include "crystal.hpp"

class monte_carlo {
public:
    double r_max;
    crystal * crystalp;
    double beta;
    monte_carlo(crystal * c) {
        crystalp = c;
        r_max = 1;
        beta = 1;
    }

    bool step_1p(particle * p) {
        vec3 candidate;
        while (true) {
            candidate = r_max * (vec3(2,2,2).random() - vec3(1,1,1));
            if (p->cell->contains(p->pos + candidate)) break;
        }
        double old_energy = p->energy();
        double new_energy = p->energy(candidate);
        double p_accept = exp(-beta*(new_energy-old_energy));
        bool accept = (rand() / (double)RAND_MAX) < std::min(p_accept, 1.);
        if (accept) {
            p->pos = crystalp->space.clip(p->pos + candidate);
        }
        return accept;
    }

    bool step_sym(particle * p1, particle * p2) {
        vec3 candidate;
        while (true) {
            candidate = r_max * (vec3(2,2,2).random() - vec3(1,1,1));
            if (p1->cell->contains(p1->pos + candidate) &&
                p2->cell->contains(p2->pos - candidate)) break;
        }
        /*if ((p1->cell->n == vec3(2, 2, 2) && p1->cell->basis == 0) ||
            (p2->cell->n == vec3(2, 2, 2) && p2->cell->basis == 0)) {
            candidate = candidate / 20;
        }*/
        double old_energy = crystalp->two_particle_energy(p1, p2);
        double new_energy = crystalp->two_particle_energy(p1, p2, candidate, -candidate);
        double p_accept = exp(-beta*(new_energy-old_energy));
        bool accept = (rand() / (double)RAND_MAX) < std::min(p_accept, 1.);
        if (accept) {
            p1->pos = crystalp->space.clip(p1->pos + candidate);
            p2->pos = crystalp->space.clip(p2->pos - candidate);
        }
        return accept;
    }

    double sweep_1p(int times=1) {
        int naccept = 0;
        for (int time = 0; time < times; time++) {
            for (particle * p : crystalp->particles) {
                naccept += step_1p(p);
            }
        }
        return (double)naccept / times / crystalp->particles.size();
    }

    double sweep_sym(int times=1) {
        int naccept = 0;
        for (int time = 0; time < times; time++) {
            int idxp1 = 0;
            for (particle * p : crystalp->particles) {
                int idxp2 = rand() % (crystalp->particles.size() - 1);
                if (idxp2 >= idxp1) idxp2 += 1;
                naccept += step_sym(p, crystalp->particles[idxp2]);
                idxp1 += 1;
            }
        }
        return (double)naccept / times / crystalp->particles.size();
    }

    double sweep(int ntimes, bool sym=true) {
        if (sym) return sweep_sym(ntimes);
        else return sweep_1p(ntimes);
    }

    void train(bool sym=true, double pacc_goal=0.3, double r_from=0.0001, double r_to=3, int width=50, int height=3, int nsweeps=5) {
        double r_best = r_from;
        double pacc_dist_best = 1;
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                double r_test = r_from + (double)i/width*(r_to-r_from);
                r_max = r_test;
                double pacc_dist = abs(sweep(nsweeps, sym) - pacc_goal);
                if (pacc_dist < pacc_dist_best) {
                    r_best = r_test;
                    pacc_dist_best = pacc_dist;
                }
            }
            double range = (r_to - r_from) / width;
            r_from = r_best - range/2;
            if (r_from < 0) r_from = 1e-8;
            r_to = r_best + range/2;
        }
        r_max = r_best;
    }
};

#endif
