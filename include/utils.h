#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <cassert>

/* Constants */
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

/* Basic Classes Claims */
class Vec;
class Interval;

/* Aliases */
using point_t = Vec;
using color_t = Vec;
using direction_t = Vec;
using tick_t = double;
using texture_t = Vec;

/* Utility Functions */
inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

// Returns a random real in [0,1).
inline double random_double() {
    return std::rand() / (RAND_MAX + 1.0);
}

// Returns a random real in [min,max).
inline double random_double(double min, double max) {
    return min + (max-min)*random_double();
}

// Convert linear color to gamma representation read by computer
inline double linear_to_gamma(double linear_component)
{
    if (linear_component > 0)
        return std::sqrt(linear_component);

    return 0;
}

/* Basic Classes Definitions */
class Vec {
  public:
    double e[3];

    Vec() : e{0,0,0} {}
    Vec(double e0, double e1, double e2) : e{e0, e1, e2} {}

    void reset() {
        e[0] = 0;
        e[1] = 0;
        e[2] = 0;
    }

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    Vec operator-() const { return Vec(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double& operator[](int i) { return e[i]; }

    Vec& operator+=(const Vec& v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    Vec& operator*=(double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    Vec& operator/=(double t) {
        assert(t != 0);
        return *this *= 1/t;
    }

    // Return the length of a vector
    double length() const {
        return std::sqrt(length_squared());
    }

    // Return the squared length of a vector
    double length_squared() const {
        return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
    }

    // Return true if the vector is close to zero in all dimensions: inside (-1e-8, 1e-8)^3
    bool near_zero() const {
        auto s = 1e-8;
        return (std::fabs(e[0]) < s) && (std::fabs(e[1]) < s) && (std::fabs(e[2]) < s);
    }
    
    // Return a random vector inside [0, 1)^3
    static Vec random() {
        return Vec(random_double(), random_double(), random_double());
    }

    // Return a random vector inside [min. max]^3
    static Vec random(double min, double max) {
        return Vec(random_double(min,max), random_double(min,max), random_double(min,max));
    }
};

inline Vec operator+(const Vec& u, const Vec& v) {
    return Vec(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline Vec operator-(const Vec& u, const Vec& v) {
    return Vec(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline Vec operator*(const Vec& u, const Vec& v) {
    return Vec(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline Vec operator*(double t, const Vec& v) {
    return Vec(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline Vec operator*(const Vec& v, double t) {
    return t * v;
}

inline Vec operator/(const Vec& v, double t) {
    assert(t != 0);
    return (1/t) * v;
}

inline double dot(const Vec& u, const Vec& v) {
    return u.e[0] * v.e[0]
         + u.e[1] * v.e[1]
         + u.e[2] * v.e[2];
}

inline Vec cross(const Vec& u, const Vec& v) {
    return Vec(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

// Return the unit vector in direction v
inline Vec unit_vector(const Vec& v) {
    assert(v.length() != 0);
    return v / v.length();
}

// Return a random vector inside a 0-centered unit disk
inline Vec random_in_unit_disk() {
    while (true) {
        auto p = Vec(random_double(-1,1), random_double(-1,1), 0);
        if (p.length_squared() < 1)
            return p;
    }
}

inline Vec random_unit_vector() {
    while (true) {
        auto p = Vec::random(-1,1);
        auto lensq = p.length_squared();
        if (1e-160 < lensq && lensq <= 1.0)
            return p / sqrt(lensq);
    }
}

inline Vec random_on_hemisphere(const Vec& normal) {
    Vec on_unit_sphere = random_unit_vector();
    if (dot(on_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return on_unit_sphere;
    else
        return -on_unit_sphere;
}

// v: in ray
// n: normal; dot(n, v) >= 0
inline Vec reflect(const Vec& v, const Vec& n) {
    return v - 2*dot(v,n)*n;
}

// uv: in ray
// n: normal; same orientation as in ray
// etai_over_etat: in/out
inline Vec refract(const Vec& uv, const Vec& n, double etai_over_etat) {
    
    auto cos_theta = std::fmin(dot(-uv, n), 1.0);
    Vec r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    Vec r_out_parallel = -std::sqrt(std::fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

class Interval {
  public:
    double min, max;

    // Default interval is empty
    Interval() : min(+infinity), max(-infinity) {} 

    Interval(double min, double max) : min(min), max(max) {}

    // Create the interval tightly enclosing the two input intervals.
    Interval(const Interval& a, const Interval& b) {
        min = a.min <= b.min ? a.min : b.min;
        max = a.max >= b.max ? a.max : b.max;
    }

    double size() const {
        return max - min;
    }

    // Whether x $\in$ Interval
    bool contains(double x) const {
        return min <= x && x <= max;
    }

    // Whether x $\in$ Interval and x is not the boundary
    bool surrounds(double x) const {
        return min < x && x < max;
    }

    // Clamp x into [min, max]
    double clamp(double x) const {
        if (x < min) return min;
        if (x > max) return max;
        return x;
    }

    // Expand current interval by `delta` and return a new expanded interval
    Interval expand(double delta) const {
        auto padding = delta/2;
        return Interval(min - padding, max + padding);
    }

    static const Interval empty, universe;
};

const Interval Interval::empty    = Interval(+infinity, -infinity);

const Interval Interval::universe = Interval(-infinity, +infinity);

// Translate the interval `ival` by `displacement`
Interval operator+(const Interval& ival, double displacement) {
    return Interval(ival.min + displacement, ival.max + displacement);
}

// Translate the interval `ival` by `displacement`
Interval operator+(double displacement, const Interval& ival) {
    return ival + displacement;
}

// Write a pixel color to a line of a file: "r g b\n"
color_t get_color(const color_t& pixel_color) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Apply a linear to gamma transform for gamma 2
    r = linear_to_gamma(r);
    g = linear_to_gamma(g);
    b = linear_to_gamma(b);

    // Translate the [0,1] component values to the byte range [0,255].
    static const Interval intensity(0.000, 0.999);
    int rbyte = int(256 * intensity.clamp(r));
    int gbyte = int(256 * intensity.clamp(g));
    int bbyte = int(256 * intensity.clamp(b));

    return color_t(rbyte, gbyte, bbyte);
}


#endif
