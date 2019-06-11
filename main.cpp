// RUN set -ex
// RUN g++ main.cpp -g -O3 -Wall
// RUN rm -fr sim
// RUN mkdir sim
// RUN time ./a.out
// RUN # sh -c 'opengl-mol sim/opengl.*'

#define NDEBUG // for a 2x speed up, disables all assert() macros

// one of these 2 options
// #define BCC_INT
#define HEXAGONAL_VAC

#include <iomanip>

#define PRINT_VAR(x) std::cout << #x" => " << (x) << std::endl

#include "crystal.hpp"
#include "monte_carlo.hpp"

#ifdef BCC_INT
crystal * crystal = crystal::bcc(4, 4, 4);
#endif

#ifdef HEXAGONAL_VAC
crystal * crystal = crystal::hexagonal(6, 6, 10);
#endif

monte_carlo monte_carlo(crystal);

void configure(double kbt_eta, double rho_sigma3) {
    double rho = crystal->density();
    crystal->potential_epsilon = 1;
    crystal->potential_sigma = pow(rho_sigma3 / rho, 1./3);
    monte_carlo.beta = 1./(kbt_eta*crystal->potential_epsilon);
}

std::string seqfn(int i) {
    std::ostringstream s;
    s << "sim/opengl." << std::setfill('0') << std::setw(4) << i;
    return s.str();
}

int main() {
    std::ofstream log_stream("sim/log");
    srand(0);
    crystal->wigner_seitz_constraint = true;
#ifdef BCC_INT
    configure(0.002, 2.5);
    lattice_cell * mid = crystal->get_cell(2, 2, 2, 0);
    particle * in = mid->interstitial(vec3(0.3, 0.3, 0.3));
    if (crystal->wigner_seitz_constraint) {
        mid->particles[0]->color = 2;
        in->color = 2;
    }
#endif
#ifdef HEXAGONAL_VAC
    configure(0.001, 4.0);
    lattice_cell * mid = crystal->get_cell(3, 3, 5, 0);
    mid->vacancy();
#endif
    monte_carlo.train();
    for (int i = 0; i < 30; i++) {
        monte_carlo.sweep_sym(100);
        crystal->write(seqfn(i));
        crystal->log(i, log_stream);
    }
    PRINT_VAR(crystal->potential_epsilon);
    PRINT_VAR(crystal->potential_sigma);
    PRINT_VAR(monte_carlo.beta);
    log_stream.close();
}
