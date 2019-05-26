#ifndef MATRIX3_HPP
#define MATRIX3_HPP

#include "vec3.hpp"

class matrix3 {
public:
    vec3 row1;
    vec3 row2;
    vec3 row3;
    vec3 col1() const { return vec3(row1.x, row2.x, row3.x); }
    vec3 col2() const { return vec3(row1.y, row2.y, row3.y); }
    vec3 col3() const { return vec3(row1.z, row2.z, row3.z); }
    matrix3 transpose() const {
        return from_cols(row1, row2, row3);
    }
    matrix3 invert() const {
        /* http://mathworld.wolfram.com/MatrixInverse.html */
        double a11 = row1.x, a12 = row1.y, a13 = row1.z;
        double a21 = row2.x, a22 = row2.y, a23 = row2.z;
        double a31 = row3.x, a32 = row3.y, a33 = row3.z;
        double det3 = det();
        auto d = [&](double a, double b, double c, double d) {
            return (a*d - b*c) / det3;
        };
        return from_rows(
                vec3(d(a22, a23, a32, a33), d(a13, a12, a33, a32), d(a12, a13, a22, a23)),
                vec3(d(a23, a21, a33, a31), d(a11, a13, a31, a33), d(a13, a11, a23, a21)),
                vec3(d(a21, a22, a31, a32), d(a12, a11, a32, a31), d(a11, a12, a21, a22))
        );
    }
    double det() const {
        /* https://en.wikipedia.org/wiki/Determinant */
        double a = row1.x, b = row1.y, c = row1.z;
        double d = row2.x, e = row2.y, f = row2.z;
        double g = row3.x, h = row3.y, i = row3.z;
        return a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;
    }
    static matrix3 from_rows(vec3 row1, vec3 row2, vec3 row3) {
        matrix3 m;
        m.row1 = row1;
        m.row2 = row2;
        m.row3 = row3;
        return m;
    }
    static matrix3 from_cols(vec3 col1, vec3 col2, vec3 col3) {
        matrix3 m;
        m.row1 = vec3(col1.x, col2.x, col3.x);
        m.row2 = vec3(col1.y, col2.y, col3.y);
        m.row3 = vec3(col1.z, col2.z, col3.z);
        return m;
    }
};

vec3 operator*(const matrix3 & a, const vec3 & b) {
    return vec3(a.row1.dot(b), a.row2.dot(b), a.row3.dot(b));
}

matrix3 operator*(const matrix3 & a, const matrix3 & b) {
    return matrix3::from_rows(
            vec3(a.row1*b.col1(), a.row1*b.col2(), a.row1*b.col3()),
            vec3(a.row2*b.col1(), a.row2*b.col2(), a.row2*b.col3()),
            vec3(a.row3*b.col1(), a.row3*b.col2(), a.row3*b.col3()));
}

std::ostream& operator<<(std::ostream &strm, const matrix3 &a) {
    return strm << "matrix3(\n  "
        << a.row1 << ",\n  "
        << a.row2 << ",\n  "
        << a.row3 << ")\n";
}

#endif
