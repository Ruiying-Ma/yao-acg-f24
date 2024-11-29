#ifndef RAY_H
#define RAY_H

#include "utils.h"

class Ray {
    public:
    Ray() {}

    Ray(const point_t& origin, const direction_t& direction) : orig_(origin), dir_(direction), time_(0.0) {}

    Ray(const point_t& origin, const direction_t& direction, tick_t time) : orig_(origin), dir_(direction), time_(time) {}

    const point_t& origin() const  { return orig_; }
    const direction_t& direction() const { return dir_; }
    tick_t time() const {return time_;}

    // Return the point on the ray which is of distance `t` from the ray origin
    point_t at(double t) const {
        return orig_ + t*dir_;
    }

  private:
    point_t orig_;
    direction_t dir_;
    tick_t time_; // a double in [0, 1)
};

#endif