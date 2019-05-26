#ifndef LATTICE_CELL_HPP
#define LATTICE_CELL_HPP

#include <cassert>
#include <vector>

#include "vec3.hpp"
#include "particle.hpp"

class particle;
class crystal;

class lattice_cell {
public:
    vec3 n /* int n1, n2, n3 */;
    vec3 nb /* periodic extent space position */;
    int basis /* int n4 */;
    vec3 center;
    std::vector<lattice_cell*> nearest_neighbours;
    std::vector<particle*> particles;
    crystal * owner;

    bool contains(const vec3 & pos) const;
    particle * interstitial(vec3 offset);
};

#endif
