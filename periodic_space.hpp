#ifndef PERIODIC_SPACE_HPP
#define PERIODIC_SPACE_HPP

#include <cassert>

#include "vec3.hpp"
#include "matrix3.hpp"

class periodic_space {
    /* this class implements coordinate transformations in a periodic lattice */
    matrix3 mat_project; /* from extent space to periodic space, columns are unit cell p1, p2 and p3 */
    matrix3 mat_project_inv; /* from periodic space to extent space */
    vec3 extent; /* real space ([0,extent.x], [0,extent.y], [0,extent.z]) */
public:
    periodic_space(matrix3 mat_project, vec3 extent) :
        mat_project(mat_project), mat_project_inv(mat_project.invert()), extent(extent) {
            assert((mat_project * (mat_project_inv * vec3(1,0,0))).close_to(vec3(1,0,0)));
            assert((mat_project * (mat_project_inv * vec3(1,-5,3))).close_to(vec3(1,-5,3)));
            assert(clip(vec3()).close_to(vec3()));
            assert(clip(mat_project.col1()/2).close_to(mat_project.col1()/2));
            assert(clip(mat_project.col2()/2).close_to(mat_project.col2()/2));
            assert(clip(mat_project.col3()/2).close_to(mat_project.col3()/2));
    }
    vec3 project(const vec3 & a) const { return mat_project * a; }
    vec3 difference(vec3 a, vec3 b) const {
        vec3 u = mat_project_inv * (b - a);
        while (u.x >= extent.x / 2) u.x -= extent.x;
        while (u.y >= extent.y / 2) u.y -= extent.y;
        while (u.z >= extent.z / 2) u.z -= extent.z;
        while (u.x < -extent.x / 2) u.x += extent.x;
        while (u.y < -extent.y / 2) u.y += extent.y;
        while (u.z < -extent.z / 2) u.z += extent.z;
        return mat_project * u;
    }
    double distance(vec3 a, vec3 b) const {
        return difference(a, b).length();
    }
    vec3 clip(const vec3 & a) const {
        vec3 u = (mat_project_inv * a).div(extent);
        u.x = u.x - floor(u.x);
        u.y = u.y - floor(u.y);
        u.z = u.z - floor(u.z);
        return (mat_project * u.mul(extent));
    }
    vec3 image(const vec3 & a, const vec3 & b) const {
        assert(0 && "do not use");
        vec3 image(0, 0, 0);
        vec3 u = mat_project_inv * (b - a);
        while (u.x >= extent.x / 2) {
            u.x -= extent.x;
            image.x -= extent.x;
        }
        while (u.y >= extent.y / 2) {
            u.y -= extent.y;
            image.y -= extent.y;
        }
        while (u.z >= extent.z / 2) {
            u.z -= extent.z;
            image.z -= extent.z;
        }
        while (u.x < -extent.x / 2) {
            u.x += extent.x;
            image.x += extent.x;
        }
        while (u.y < -extent.y / 2) {
            u.y += extent.y;
            image.y += extent.y;
        }
        while (u.z < -extent.z / 2) {
            u.z += extent.z;
            image.z += extent.z;
        }
        return mat_project * image;
    }

    vec3 p1() const { return extent.x * mat_project.col1(); }
    vec3 p2() const { return extent.y * mat_project.col2(); }
    vec3 p3() const { return extent.z * mat_project.col3(); }
    double volume() const {
        return abs(p1() * (p2().cross(p3())));
    }
};

#endif
