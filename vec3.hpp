#ifndef VEC3_HPP
#define VEC3_HPP

#include <math.h>
#include <iostream>

class vec3 {
public:
    double x = 0;
    double y = 0;
    double z = 0;
    vec3() {}
    vec3(double x, double y, double z) : x(x), y(y), z(z) { }
    vec3 operator+(const vec3 & other) const { return vec3(x + other.x, y + other.y, z + other.z); }
    vec3 operator-(const vec3 & other) const { return vec3(x - other.x, y - other.y, z - other.z); }
    double operator*(const vec3 & other) const { return x * other.x + y * other.y + z * other.z; }
    vec3 operator*(double a) const { return vec3(x * a, y * a, z * a); }
    double dot(const vec3 & other) const { return *this * other; }
    vec3 operator/(double a) const { return vec3(x / a, y / a, z / a); }
    vec3 operator-() const { return vec3(-x, -y, -z); }
    vec3 operator+() const { return *this; }
    vec3 & operator+=(const vec3 & other) { x += other.x; y += other.y; z += other.z; return *this; }
    vec3 & operator-=(const vec3 & other) { x -= other.x; y -= other.y; z -= other.z; return *this; }
    double length() const { return sqrt(x*x + y*y + z*z); }
    vec3 unit() const { return *this / length(); }
    double dist(const vec3 & to) const { return (to - *this).length(); }
    vec3 octant(const vec3 & point) const { return vec3(point.x < x ? -1 : 1, point.y < y ? -1 : 1, point.z < z ? -1 : 1); }
    vec3 mul(const vec3 & other) const { return vec3(x*other.x, y*other.y, z*other.z); }
    vec3 div(const vec3 & other) const { return vec3(x/other.x, y/other.y, z/other.z); }
    vec3 mul(double xx, double yy, double zz) const { return vec3(x*xx, y*yy, z*zz); }
    double cos(const vec3 & other) const { return dot(other) / length() / other.length(); }
    double acos(const vec3 & other) const { return std::acos(cos(other)); }
    vec3 cross(const vec3 & o) const { return vec3(y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x); }
    vec3 random() const { return (vec3(rand(), rand(), rand()) / RAND_MAX).mul(*this); }
    bool operator==(const vec3 & other) const { return x == other.x && y == other.y && z == other.z; }
    bool close_to(const vec3 & other, double tol=1e-4) const {
        return abs(x - other.x) < tol && abs(y - other.y) < tol && abs(z - other.z) < tol;
    }
};

inline vec3 operator*(double a, const vec3 & self) { return self * a; }
inline vec3 operator/(double a, const vec3 & self) { return self / a; }

std::ostream& operator<<(std::ostream &strm, const vec3 &a) {
    return strm << "vec3(" << a.x << "," << a.y << "," << a.z << ")";
}

#endif
