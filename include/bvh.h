#ifndef BVH_H
#define BVH_H

#include "hittable.h"
#include <algorithm>


class BVHNode : public Hittable {
  public:
    BVHNode(HittableVector list) : BVHNode(list.hittable_vec, 0, list.hittable_vec.size()) {
        // There's a C++ subtlety here. This constructor (without span indices) creates an
        // implicit copy of the Hittable list, which we will modify. The lifetime of the copied
        // list only extends until this constructor exits. That's OK, because we only need to
        // persist the resulting bounding volume hierarchy.
    }

    // Build the bounding box of the span [start, end) of source objects.
    BVHNode(std::vector<std::shared_ptr<Hittable>>& hittable_vec, size_t start, size_t end) {
        bbox_ = BoundingBox::empty;
        for (size_t hittable_index=start; hittable_index < end; hittable_index++)
            bbox_ = BoundingBox(bbox_, hittable_vec[hittable_index]->bounding_box());

        int axis = bbox_.longest_axis();

        auto comparator = (axis == 0) ? box_x_compare
                        : (axis == 1) ? box_y_compare
                                      : box_z_compare;

        size_t span = end - start;

        if (span == 1) {
            left_hittable_ = right_hittable_ = hittable_vec[start];
        } else if (span == 2) {
            left_hittable_ = hittable_vec[start];
            right_hittable_ = hittable_vec[start+1];
        } else {
            std::sort(std::begin(hittable_vec) + start, std::begin(hittable_vec) + end, comparator);

            auto mid = start + span/2;
            left_hittable_ = std::make_shared<BVHNode>(hittable_vec, start, mid);
            right_hittable_ = std::make_shared<BVHNode>(hittable_vec, mid, end);
        }
    }

    bool hit(const Ray& r, Interval dist_range, HitRecord& hit_rec) const override {
        if (!bbox_.hit(r, dist_range))
            return false;

        bool hit_left = left_hittable_->hit(r, dist_range, hit_rec);
        bool hit_right = right_hittable_->hit(r, Interval(dist_range.min, hit_left ? hit_rec.dist : dist_range.max), hit_rec);

        return hit_left || hit_right;
    }

    BoundingBox bounding_box() const override { return bbox_; }

  private:
    std::shared_ptr<Hittable> left_hittable_;
    std::shared_ptr<Hittable> right_hittable_;
    BoundingBox bbox_;

    // Compare the min margin of the bounding boxes of two hittables along the given axis
    static bool box_compare(
        const std::shared_ptr<Hittable> a, const std::shared_ptr<Hittable> b, int axis_index
    ) {
        auto a_axis_interval = a->bounding_box().axis_interval(axis_index);
        auto b_axis_interval = b->bounding_box().axis_interval(axis_index);
        return a_axis_interval.min < b_axis_interval.min;
    }

    static bool box_x_compare (const std::shared_ptr<Hittable> a, const std::shared_ptr<Hittable> b) {
        return box_compare(a, b, 0);
    }

    static bool box_y_compare (const std::shared_ptr<Hittable> a, const std::shared_ptr<Hittable> b) {
        return box_compare(a, b, 1);
    }

    static bool box_z_compare (const std::shared_ptr<Hittable> a, const std::shared_ptr<Hittable> b) {
        return box_compare(a, b, 2);
    }
};

#endif