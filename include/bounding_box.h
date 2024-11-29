#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include "ray.h"

// A 3d bounding box used for BVH acceleration
class BoundingBox {
    public:
    Interval x, y, z;

    // The default BoundingBox is empty, since intervals are empty by default.
    BoundingBox() {} 

    BoundingBox(const Interval& x, const Interval& y, const Interval& z)
      : x(x), y(y), z(z)
    {
        pad_to_minimums_();
    }

    // Treat the two points a and b as extrema (top_left & bottom_right) for the bounding box, so we don't require a particular minimum/maximum coordinate order.
    BoundingBox(const point_t& a, const point_t& b) {
        x = (a[0] <= b[0]) ? Interval(a[0], b[0]) : Interval(b[0], a[0]);
        y = (a[1] <= b[1]) ? Interval(a[1], b[1]) : Interval(b[1], a[1]);
        z = (a[2] <= b[2]) ? Interval(a[2], b[2]) : Interval(b[2], a[2]);

        pad_to_minimums_();
    }

    // Create the bounding box that tightly enclosing the two bounding boxes
    BoundingBox(const BoundingBox& box0, const BoundingBox& box1) {
        x = Interval(box0.x, box1.x);
        y = Interval(box0.y, box1.y);
        z = Interval(box0.z, box1.z);
    }

    // Return the range of this bounding box in direction `n` (n = 0, 1, 2)
    // Default as n=0 (x)
    const Interval& axis_interval(int n) const {
        if (n == 1) return y;
        if (n == 2) return z;
        return x;
    }

    // Return whether the ray intersects this bounding box
    // hit_range: the intersection of the the hit ranges of the three directions
    bool hit(const Ray& r, Interval hit_range) const {
        const point_t& ray_orig = r.origin();
        const direction_t& ray_dir  = r.direction();

        for (int axis = 0; axis < 3; axis++) {
            const Interval& ax = axis_interval(axis);
            const double adinv = 1.0 / ray_dir[axis];

            auto t0 = (ax.min - ray_orig[axis]) * adinv;
            auto t1 = (ax.max - ray_orig[axis]) * adinv;

            if (t0 < t1) {
                if (t0 > hit_range.min) hit_range.min = t0;
                if (t1 < hit_range.max) hit_range.max = t1;
            } else {
                if (t1 > hit_range.min) hit_range.min = t1;
                if (t0 < hit_range.max) hit_range.max = t0;
            }

            if (hit_range.max <= hit_range.min)
                return false;
        }
        return true;
    }

    // Returns the index of the longest interval of the bounding box.
    int longest_axis() const {
        if (x.size() > y.size())
            return x.size() > z.size() ? 0 : 2;
        else
            return y.size() > z.size() ? 1 : 2;
    }

    static const BoundingBox empty, universe;

  private:
    // Adjust the BoundingBox so that no side is narrower than some delta, padding if necessary.
    void pad_to_minimums_() {
        double delta = 0.0001;
        if (x.size() < delta) x = x.expand(delta);
        if (y.size() < delta) y = y.expand(delta);
        if (z.size() < delta) z = z.expand(delta);
    }
};

const BoundingBox BoundingBox::empty = BoundingBox(Interval::empty,    Interval::empty, Interval::empty);

const BoundingBox BoundingBox::universe = BoundingBox(Interval::universe, Interval::universe, Interval::universe);

// Translate the bounding box `bbox` by vector `offset`
BoundingBox operator+(const BoundingBox& bbox, const Vec& offset) {
    return BoundingBox(bbox.x + offset.x(), bbox.y + offset.y(), bbox.z + offset.z());
}

// Translate the bounding box `bbox` by vector `offset`
BoundingBox operator+(const Vec& offset, const BoundingBox& bbox) {
    return bbox + offset;
}

#endif