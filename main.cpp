// RUN set -ex
// RUN g++ main.cpp -g -O3
// RUN rm -fr sim
// RUN mkdir sim
// RUN time ./a.out
// RUN sh -c 'opengl-mol sim/*'

#include <iomanip>

#define PRINT_VAR(x) std::cout << #x" = " << (x) << std::endl

#include "crystal.hpp"
#include "monte_carlo.hpp"

crystal * crystal = crystal::bcc();
monte_carlo monte_carlo(crystal);

void configure(double kbt_eta, double rho_sigma3) {
    double rho = crystal->density();
    crystal->potential_epsilon = 1;
    crystal->potential_sigma = pow(rho_sigma3 / rho, 1./3);
    monte_carlo.beta = 1./(kbt_eta*crystal->potential_epsilon);
}

std::string seqfn(int i){ 
    std::ostringstream s;
    s << "sim/opengl."
      << std::setfill('0') << std::setw(4)
      << i;
    return s.str();
}

int main() {
    srand(0);
    configure(0.003, 2.5);
    crystal->wigner_seitz_constraint = false;
    monte_carlo.train();
    lattice_cell * mid = crystal->get_cell(2, 2, 2, 0);
    particle * in = mid->interstitial(vec3(0.3, 0.3, 0.3));
    if (crystal->wigner_seitz_constraint) {
        mid->particles[0]->color = 2;
        in->color = 2;
    }
    monte_carlo.train();
    for (int i = 0; i < 100; i++) {
        monte_carlo.sweep_sym(100);
        crystal->write(seqfn(i));
    }
    PRINT_VAR(crystal->potential_epsilon);
    PRINT_VAR(crystal->potential_sigma);
    PRINT_VAR(monte_carlo.beta);
}
