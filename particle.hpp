#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "vec3.hpp"

class lattice_cell;
class crystal;

class particle {
public:
    vec3 pos;
    lattice_cell * cell;
    crystal * owner;

    double energy(vec3 shift=vec3(0, 0, 0));
};

#endif
