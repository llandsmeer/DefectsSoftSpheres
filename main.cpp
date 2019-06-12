// RUN set -ex
// RUN g++ main.cpp -g -O3 -Wall
// RUN rm -fr sim
// RUN mkdir sim
// RUN time ./a.out
// RUN # sh -c 'opengl-mol sim/opengl.*'

// #define NDEBUG // for a 2x speed up, disables all assert() macros

// one of these 2 options
#define HERTZ_BCC_INT
// #define HERTZ_HEX_VAC

#include <math.h>
#include <iomanip>

#define PRINT_VAR(x) std::cout << #x" => " << (x) << std::endl

#include "crystal.hpp"
#include "monte_carlo.hpp"
#include "axis_offsets.hpp"
#include "bcc_offsets.hpp"

#ifdef HERTZ_BCC_INT
crystal * crystal = crystal::bcc(7, 7, 7);
#endif

#ifdef HERTZ_HEX_VAC
crystal * crystal = crystal::hexagonal(6, 6, 20);
axis_offsets axis_offsets(crystal, vec3(0, 0, 1));
#endif

monte_carlo monte_carlo(crystal);

void configure_hertz(double kbt_eta, double rho_sigma3) {
    crystal->potential_epsilon = 1;
    crystal->potential_sigma = pow(rho_sigma3 / crystal->density(), 1./3);
    monte_carlo.beta = 1./(kbt_eta*crystal->potential_epsilon);
}

void configure_star(double packing_fraction, double one_over_f) {
    crystal->potential_type = crystal::potential_type::STAR;
    crystal->potential_sigma = pow(M_PI / (6. * crystal->density() * packing_fraction), 1/3.);
    crystal->potential_epsilon = 1. / one_over_f;
    monte_carlo.beta = 1;
}

std::string seqfn(int i) {
    std::ostringstream s;
    s << "sim/opengl." << std::setfill('0') << std::setw(4) << i;
    return s.str();
}

int main() {
#ifdef NDEBUG
    std::cout << "WARNING NDEBUG IS DEFINED\n";
#endif
    srand(0);
    crystal->wigner_seitz_constraint = true;
#ifdef HERTZ_BCC_INT
    configure_hertz(0.002, 2.5);
    std::ofstream log_stream("sim/bcc_offsets");
    lattice_cell * mid = crystal->get_cell(4, 4, 4, 0);
    particle * in = mid->interstitial(vec3(0.3, 0.3, 0.3));
    if (crystal->wigner_seitz_constraint) {
        mid->particles[0]->color = 2;
        in->color = 2;
    }
    bcc_offsets axis_offsets(crystal, mid);
#endif
#ifdef HERTZ_HEX_VAC
    configure_hertz(0.001, 4.0);
    lattice_cell * mid = crystal->get_cell(3, 3, 5, 0);
    mid->vacancy();
    std::ofstream log_stream("sim/hex_offsets");
    for (int i = 0; i < 20; i++) {
        lattice_cell * lc = crystal->get_cell(3, 3, i, 0);
        for (particle * p : lc->particles) {
            axis_offsets.add_particle(p);
        }
    }
#endif
    monte_carlo.train();
    for (int i = 0; i < 100; i++) {
        monte_carlo.sweep_sym(100);
        crystal->write(seqfn(i));
        axis_offsets.measure();
        axis_offsets.write(log_stream);
    }
    PRINT_VAR(crystal->potential_epsilon);
    PRINT_VAR(crystal->potential_sigma);
    PRINT_VAR(monte_carlo.beta);
    log_stream.close();
}
