#ifndef LATTICE_DEFINITION_HPP
#define LATTICE_DEFINITION_HPP

#include <math.h>
#include <vector>

#include "vec3.hpp"

class lattice_definition {
private:
    double cos_alpha_prime() const {
        return
            (cos(alpha) - cos(beta) * cos(gamma)) / 
            (sin(beta)*sin(gamma));
    }

    double sin_alpha_prime() const {
        double c2a = cos(alpha) * cos(alpha);
        double c2b = cos(beta) * cos(beta);
        double s2g = sin(gamma) * sin(gamma);
        return
            sqrt(-c2a-c2b+s2g+2*cos(alpha)*cos(beta)*cos(gamma)) /
            (sin(beta) * sin(gamma));
    };
public:
    double a;
    double b;
    double c;
    double alpha;
    double beta;
    double gamma;
    std::vector<vec3> basis_vectors;

    vec3 p1() const {
        return vec3(a, 0, 0);
    }

    vec3 p2() const {
        double p2x = b*cos(gamma);
        double p2y = b*sin(gamma);
        return vec3(
                abs(p2x) < 1e-10 ? 0 : p2x,
                abs(p2y) < 1e-10 ? 0 : p2y,
                0);
    }

    vec3 p3() const {
        double p3x = c*cos(beta);
        double p3y = c*cos_alpha_prime()*sin(beta);
        double p3z = c*sin_alpha_prime()*sin(beta);
        return vec3(
                abs(p3x) < 1e-10? 0 : p3x,
                abs(p3y) < 1e-10? 0 : p3y,
                abs(p3z) < 1e-10? 0 : p3z);
    }

    double volume() const {
        return a * b * c * sin_alpha_prime()*sin(beta)*sin(gamma);
    }

    static lattice_definition simple_cubic(double a=1) {
        lattice_definition lattice;
        lattice.a = a;
        lattice.b = a;
        lattice.c = a;
        lattice.alpha = M_PI / 2;
        lattice.beta = M_PI / 2;
        lattice.gamma = M_PI / 2;
        lattice.basis_vectors.push_back(vec3());
        return lattice;
    }

    static lattice_definition face_centered_cubic(double a=1) {
        lattice_definition lattice = simple_cubic(a);
        lattice.basis_vectors.push_back(vec3(0.5, 0.5, 0.0));
        lattice.basis_vectors.push_back(vec3(0.0, 0.5, 0.5));
        lattice.basis_vectors.push_back(vec3(0.5, 0.0, 0.5));
        return lattice;
    }

    static lattice_definition body_centered_cubic(double a=1) {
        lattice_definition lattice = simple_cubic(a);
        lattice.basis_vectors.push_back(vec3(0.5, 0.5, 0.5));
        return lattice;
    }

    static lattice_definition body_centered_tetragonal(double a, double c_over_a) {
        lattice_definition lattice;
        lattice.a = a;
        lattice.b = a;
        lattice.c = a * c_over_a;
        lattice.alpha = acos(1 / (2*c_over_a));
        lattice.beta = lattice.alpha;
        lattice.gamma = M_PI / 2;
        lattice.basis_vectors.push_back(vec3());
        return lattice;
    }

    static lattice_definition hexagonal(double a, double c_over_a) {
        lattice_definition lattice;
        lattice.a = a;
        lattice.b = a;
        lattice.c = a * c_over_a;
        lattice.alpha = M_PI / 2;
        lattice.beta = M_PI / 2;
        lattice.gamma = M_PI / 3;
        lattice.basis_vectors.push_back(vec3());
        return lattice;
    }

    static lattice_definition body_centered_orthorhombic(double a, double b_over_a, double c_over_a) {
        lattice_definition lattice;
        lattice.a = a;
        lattice.b = a * b_over_a;
        lattice.c = a * c_over_a;
        lattice.alpha = M_PI / 2;
        lattice.beta = M_PI / 2;
        lattice.gamma = M_PI / 2;
        lattice.basis_vectors.push_back(vec3());
        lattice.basis_vectors.push_back(vec3(0.5, 0.5, 0.5));
        return lattice;
    }

    static lattice_definition diamond(double a) {
        /* https://en.wikipedia.org/wiki/Diamond_cubic#Mathematical_structure */
        lattice_definition lattice = simple_cubic(a);
        lattice.basis_vectors.push_back(vec3(0, 2, 2)/4);
        lattice.basis_vectors.push_back(vec3(2, 0, 2)/4);
        lattice.basis_vectors.push_back(vec3(2, 2, 0)/4);
        lattice.basis_vectors.push_back(vec3(3, 3, 3)/4);
        lattice.basis_vectors.push_back(vec3(3, 1, 1)/4);
        lattice.basis_vectors.push_back(vec3(1, 3, 1)/4);
        lattice.basis_vectors.push_back(vec3(1, 1, 3)/4);
        return lattice;
    }
};

#endif
