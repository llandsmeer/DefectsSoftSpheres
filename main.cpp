// RUN set -e
// RUN echo ===== "$1" =====
// RUN if [ "$1" = test ]; then
// RUN     rm -fr sim/test
// RUN fi
// RUN mkdir sim/"$1"
// RUN cp main.cpp sim/"$1"/main.cpp
// RUN g++ main.cpp -g -O3 -Wall -std=c++17 -o sim/"$1"/main -lpthread
// RUN time sim/"$1"/main "$1"
// RUN sh -c "opengl-mol sim/$1/"'opengl.*'

#define NDEBUG // for a 2x speed up, disables all assert() macros

// one of these 2 options
// #define HERTZ_BCC_INT
#define HERTZ_FCC_INT
// #define HERTZ_SC_VAC
// #define HERTZ_HEX_VAC
// #define STAR_TEST

#include <thread>
#include <math.h>
#include <iomanip>

#define PRINT_VAR(x) std::cout << #x" => " << (x) << std::endl

#include "crystal.hpp"
#include "monte_carlo.hpp"
#include "axis_offsets.hpp"
#include "bcc_offsets.hpp"
#include "sc_offsets.hpp"

#ifdef HERTZ_SC_VAC
auto unitcell = lattice_definition::simple_cubic(3);
crystal * crystal = crystal::build(unitcell, 10);
#endif

#ifdef HERTZ_BCC_INT
auto unitcell = lattice_definition::body_centered_cubic(3);
crystal * crystal = crystal::build(unitcell, 7, 7, 7);
#endif

#ifdef HERTZ_FCC_INT
auto unitcell = lattice_definition::face_centered_cubic(3);
crystal * crystal = crystal::build(unitcell, 8);
#endif

#ifdef HERTZ_HEX_VAC
auto unitcell = lattice_definition::hexagonal(3, 0.84);
crystal * crystal = crystal::build(unitcell, 6, 6, 20);
axis_offsets axis_offsets(crystal, vec3(0, 0, 1));
#endif

#ifdef STAR_TEST
auto unitcell_bcc = lattice_definition::body_centered_cubic(3);
auto unitcell_bco = lattice_definition::body_centered_orthorhombic(1, 3.14, 1.81);
auto unitcell_diam = lattice_definition::diamond(5);
crystal * crystal = crystal::build(unitcell_bcc, 4, 4, 4, 10);
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
    s << "opengl." << std::setfill('0') << std::setw(4) << i;
    return s.str();
}

template<typename Str>
std::string path_join(Str a) {
    std::ostringstream s; s << a; return s.str(); }
template<typename StrHead, typename... StrTail>
std::string path_join(StrHead head, StrTail... tail) {
     std::ostringstream s; s << head << "/" << path_join(tail...); return s.str(); }

int main(int argc, char ** argv) {
    volatile int progress = 0;
    std::thread([&]() { for(;;) { std::this_thread::sleep_for(std::chrono::seconds(2));
        std::cerr << progress << " " << std::flush;
    }}).detach();
    if (argc != 2) {
        std::cerr << "usage: sh c [output directory in sim/]" << std::endl;
        exit(1);
    }
    std::string root = path_join("sim", argv[1]);
#ifdef NDEBUG
    std::cout << "WARNING NDEBUG IS DEFINED\n";
#endif
    srand(0);
    crystal->wigner_seitz_constraint = true;
#ifdef HERTZ_SC_VAC
    configure_hertz(0.001, 5.2);
    std::ofstream log_stream(path_join(root, "sc_offsets"));
    lattice_cell * mid = crystal->get_cell(5, 5, 5, 0);
    mid->vacancy();
    sc_offsets axis_offsets(crystal, mid);
#endif
#ifdef HERTZ_BCC_INT
    configure_hertz(0.002, 2.5);
    std::ofstream log_stream(path_join(root, "bcc_offsets"));
    lattice_cell * mid = crystal->get_cell(4, 4, 4, 0);
    particle * in = mid->interstitial(vec3(0.3, 0.3, 0.3));
    if (crystal->wigner_seitz_constraint) {
        mid->particles[0]->color = 2;
        in->color = 2;
    }
    bcc_offsets axis_offsets(crystal, mid);
#endif
#ifdef HERTZ_FCC_INT
    crystal->wigner_seitz_constraint = true;
    configure_hertz(0.002, 1.8);
    std::ofstream log_stream(path_join(root, "fcc_offsets"));
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
    std::ofstream log_stream(path_join(root, "hex_offsets"));
    for (int i = 0; i < 20; i++) {
        lattice_cell * lc = crystal->get_cell(3, 3, i, 0);
        for (particle * p : lc->particles) {
            axis_offsets.add_particle(p);
        }
    }
#endif
#ifdef STAR_TEST
    configure_star(1.25, 0.01);
    crystal->wigner_seitz_constraint = false;
    lattice_cell * mid = crystal->get_cell(0, 0, 0, 0);
    bcc_offsets axis_offsets(crystal, mid);
    std::ofstream log_stream(path_join(root, "bcc_offsets"));
#endif
    PRINT_VAR(crystal->particles.size());
    crystal->write(path_join(root, seqfn(0)));
    monte_carlo.train();
    for (int i = 0; i < 100; i++) {
        monte_carlo.sweep_sym(100);
        crystal->write(path_join(root, seqfn(i+1)));
        axis_offsets.measure();
        axis_offsets.write(log_stream);
        progress = i;
    }
    log_stream.close();
}
