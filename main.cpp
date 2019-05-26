#define PRINT_VAR(x) std::cout << #x" = " << (x) << std::endl

#include "crystal.hpp"
#include "monte_carlo.hpp"

crystal * crystal = crystal::bcc();
monte_carlo monte_carlo(crystal);

void configure(double kbt_eta, double rho_sigma3) {
    double rho = crystal->density();
    crystal->potential_epsilon = 1;
    crystal->potential_sigma = pow(rho_sigma3 / rho, 1./3);
    monte_carlo.beta = 1./kbt_eta/crystal->potential_epsilon;
}

int main() {
    srand(0);
    configure(0.002, 2.5);
    monte_carlo.train();
    lattice_cell * mid = crystal->get_cell(2, 2, 2, 0);
    crystal->write("a.opengl");
    particle * in = mid->interstitial(vec3(0.3, 0.3, 0.3));
    crystal->write("b.opengl");
    monte_carlo.train();
    PRINT_VAR(crystal->potential_epsilon);
    PRINT_VAR(crystal->potential_sigma);
    PRINT_VAR(monte_carlo.beta);
}
