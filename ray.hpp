#ifndef RAY_HPP
#define RAY_HPP

#include "vec3.hpp"

class ray {
   public:
    point3 orig;
    vec3 dir;
    double tm;

    ray() {}
    ray(const point3& origin, const vec3& direction, double time=0.0) : orig{origin}, dir{direction}, tm{time} {}

    point3 origin() const { return orig; }
    vec3 direction() const { return dir; }
    double time() const { return tm; }

    point3 at(double t) const {
        return orig + t * dir;
    }
};

#endif