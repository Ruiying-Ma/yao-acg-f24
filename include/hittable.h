#ifndef HITTABLE_H
#define HITTABLE_H

#include "material.h"
#include "bounding_box.h"
#include <vector>


class Hittable {
    public:
    virtual ~Hittable() = default;

    // Return whether the ray will hit this object
    // The distance from the ray origin to the intersection point must be inside `dist_range`
    // record: **store** the intersection record to `hit_rec`
    virtual bool hit(const Ray& r, Interval dist_range, HitRecord& hit_rec) const = 0;

    // Return the bounding box of this object
    virtual BoundingBox bounding_box() const = 0;
};

class HittableVector : public Hittable {
  public:
    std::vector<std::shared_ptr<Hittable>> hittable_vec;

    HittableVector() {}
    HittableVector(std::shared_ptr<Hittable> hittable) { push_back(hittable); }

    void clear() { hittable_vec.clear(); }

    void push_back(std::shared_ptr<Hittable> hittable) {
        hittable_vec.push_back(hittable);
        bbox_ = BoundingBox(bbox_, hittable->bounding_box());
    }

    bool hit(const Ray& r, Interval dist_range, HitRecord& hit_rec) const override {
        HitRecord cur_hit_rec;
        bool has_hit = false;
        double shortest_dist = dist_range.max;

        for (const auto& object : hittable_vec) {
            if (object->hit(r, Interval(dist_range.min, shortest_dist), cur_hit_rec)) {
                has_hit = true;
                assert(cur_hit_rec.dist <= shortest_dist);
                shortest_dist = cur_hit_rec.dist;
                hit_rec = cur_hit_rec;
            }
        }

        return has_hit;
    }

    BoundingBox bounding_box() const override { return bbox_; }

  private:
    BoundingBox bbox_;
};

class TranslateHittable : public Hittable {
  public:
    TranslateHittable(std::shared_ptr<Hittable> obj, const Vec& offset): hittable_(obj), offset_(offset)
    {
        bbox_ = obj->bounding_box() + offset;
    }

    bool hit(const Ray& r, Interval dist_range, HitRecord& hit_rec) const override {
        // Move the ray backwards by the offset
        Ray offset_r(r.origin() - offset_, r.direction(), r.time());

        // Determine whether an intersection exists along the offset ray (and if so, where)
        if (!hittable_->hit(offset_r, dist_range, hit_rec))
            return false;

        // Move the intersection point forwards by the offset
        hit_rec.ip += offset_;

        return true;
    }

    BoundingBox bounding_box() const override { return bbox_; }

  private:
    std::shared_ptr<Hittable> hittable_;
    Vec offset_;
    BoundingBox bbox_;
};

class RotateYHittable : public Hittable {
  public:
    RotateYHittable(std::shared_ptr<Hittable> obj, double angle) : hittable_(obj) {
        auto radians = degrees_to_radians(angle);
        sin_theta_ = std::sin(radians);
        cos_theta_ = std::cos(radians);
        bbox_ = obj->bounding_box();

        // Rotate bbox_ along y-axis
        point_t min( infinity,  infinity,  infinity);
        point_t max(-infinity, -infinity, -infinity);

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    auto x = i*bbox_.x.max + (1-i)*bbox_.x.min;
                    auto y = j*bbox_.y.max + (1-j)*bbox_.y.min;
                    auto z = k*bbox_.z.max + (1-k)*bbox_.z.min;

                    auto newx =  cos_theta_*x + sin_theta_*z;
                    auto newz = -sin_theta_*x + cos_theta_*z;

                    point_t tester(newx, y, newz);

                    for (int c = 0; c < 3; c++) {
                        min[c] = std::fmin(min[c], tester[c]);
                        max[c] = std::fmax(max[c], tester[c]);
                    }
                }
            }
        }

        bbox_ = BoundingBox(min, max);
    }

    bool hit(const Ray& r, Interval dist_range, HitRecord& hit_rec) const override {
        // Transform the ray from world space to object space.
        Ray rotated_r(world2local_(r.origin()), world2local_(r.direction()), r.time());

        // Determine whether an intersection exists in object space (and if so, where).
        if (!hittable_->hit(rotated_r, dist_range, hit_rec))
            return false;

        // Transform the intersection from object space back to world space.
        hit_rec.ip = local2world_(hit_rec.ip);
        hit_rec.normal = local2world_(hit_rec.normal);

        return true;
    }

    BoundingBox bounding_box() const override { return bbox_; }

  private:
    std::shared_ptr<Hittable> hittable_;
    double sin_theta_;
    double cos_theta_;
    BoundingBox bbox_;

    Vec world2local_(const Vec& v) const {
        return Vec(
            (cos_theta_ * v.x()) - (sin_theta_ * v.z()),
            v.y(),
            (sin_theta_ * v.x()) + (cos_theta_ * v.z())
        );
    }

    Vec local2world_(Vec v) const {
        return Vec(
            (cos_theta_ * v.x()) + (sin_theta_ * v.z()),
            v.y(),
            (-sin_theta_ * v.x()) + (cos_theta_ * v.z())
        );
    }
};

class Sphere : public Hittable {
    public:
    // Stationary Sphere
    Sphere(const point_t& static_center, double radius, std::shared_ptr<Material> material)
      : center_(static_center, point_t(0,0,0)), radius_(std::fmax(0,radius)), material_(material)
    {
        auto rvec = Vec(radius, radius, radius);
        bbox_ = BoundingBox(static_center - rvec, static_center + rvec);
    }

    // Moving Sphere
    Sphere(const point_t& center_from, const point_t& center_to, double radius,
           std::shared_ptr<Material> material)
      : center_(center_from, center_to - center_from), radius_(std::fmax(0,radius)), material_(material)
    {
        auto rvec = Vec(radius, radius, radius);
        BoundingBox box_from(center_.at(0) - rvec, center_.at(0) + rvec);
        BoundingBox box_to(center_.at(1) - rvec, center_.at(1) + rvec);
        bbox_ = BoundingBox(box_from, box_to);
    }
    
    bool hit(const Ray& r, Interval dist_range, HitRecord& hit_rec) const override {
        point_t current_center = center_.at(r.time());
        direction_t ray2sphere = current_center - r.origin();
        double a = r.direction().length_squared();
        double h = dot(r.direction(), ray2sphere);
        double c = ray2sphere.length_squared() - radius_*radius_;

        auto discriminant = h*h - a*c;
        if (discriminant < 0)
            return false;

        auto sqrtd = std::sqrt(discriminant);

        // Find the nearest root that lies in the acceptable range.
        auto root = (h - sqrtd) / a;
        if (!dist_range.surrounds(root)) {
            root = (h + sqrtd) / a;
            if (!dist_range.surrounds(root))
                return false;
        }

        hit_rec.dist = root;
        hit_rec.ip = r.at(hit_rec.dist);
        direction_t outward_normal = (hit_rec.ip - current_center) / radius_; // normalized
        hit_rec.set_face_normal(r, outward_normal);
        get_uv_(outward_normal, hit_rec.u, hit_rec.v);
        hit_rec.mat = material_;

        return true;
    }

    BoundingBox bounding_box() const override {
        return bbox_;
    }

    private:
    Ray center_; // time = 0
    double radius_;
    std::shared_ptr<Material> material_;
    BoundingBox bbox_;

    // Store the (u, v) of texture space, given p of world space
    // NOTE: `p` centered at the origin of the object
    // NOTE: `p` is unit vector
    static void get_uv_(const point_t& p, double& u, double& v) {
        auto theta = std::acos(-p.y());
        auto phi = std::atan2(-p.z(), p.x()) + pi;

        u = phi / (2*pi);
        v = theta / pi;
    }
};

class Quad : public Hittable {
    public:
    Quad(const point_t& orig, const Vec& off_u, const Vec& off_v, std::shared_ptr<Material> material)
      : orig_(orig), off_u_(off_u), off_v_(off_v), material_(material)
    {
        auto n = cross(off_u, off_v);
        normal_ = unit_vector(n);
        dot_orig_normal_ = dot(normal_, orig);
        repri_normal_ = n / dot(n,n);

        bbox_ = BoundingBox(BoundingBox(orig, orig + off_u + off_v), BoundingBox(orig + off_u, orig + off_v));
    }

    bool hit(const Ray& r, Interval dist_range, HitRecord& hit_rec) const override {
        double denom = dot(normal_, r.direction());

        // No hit if the ray is parallel to the plane.
        if (std::fabs(denom) < 1e-8)
            return false;

        // Return false if the hit point parameter t is outside the ray interval.
        double dist = (dot_orig_normal_ - dot(normal_, r.origin())) / denom;
        if (!dist_range.contains(dist))
            return false;

        // Determine if the hit point lies within the planar shape using its plane coordinates.
        point_t intersection = r.at(dist);
        point_t intersection_off = intersection - orig_;
        double alpha = dot(repri_normal_, cross(intersection_off, off_v_));
        auto beta = dot(repri_normal_, cross(off_u_, intersection_off));

        // Determine if the intersection point is in the finite plane
        Interval unit_interval = Interval(0, 1);
        if (!unit_interval.contains(alpha) || !unit_interval.contains(beta)) {
            return false;
        }
        hit_rec.u = alpha;
        hit_rec.v = beta;

        // Ray hits the 2D shape; set the rest of the hit record and return true.
        hit_rec.dist = dist;
        hit_rec.ip = intersection;
        hit_rec.mat = material_;
        hit_rec.set_face_normal(r, normal_);

        return true;
    }

    BoundingBox bounding_box() const override {
        return bbox_;
    }

    private:
    point_t orig_;
    Vec off_u_, off_v_;
    Vec repri_normal_; // normal_ / length(normal_)
    std::shared_ptr<Material> material_;
    BoundingBox bbox_;
    direction_t normal_; // unit length; direction follows right-hand rule
    double dot_orig_normal_;
};


class Box : public HittableVector {
    public:
    Box(const point_t& a, const point_t& b, std::shared_ptr<Material> material) {
        point_t min = point_t(std::fmin(a.x(),b.x()), std::fmin(a.y(),b.y()), std::fmin(a.z(),b.z()));
        point_t max = point_t(std::fmax(a.x(),b.x()), std::fmax(a.y(),b.y()), std::fmax(a.z(),b.z()));

        point_t dx = point_t(max.x() - min.x(), 0, 0);
        point_t dy = point_t(0, max.y() - min.y(), 0);
        point_t dz = point_t(0, 0, max.z() - min.z());

        push_back(std::make_shared<Quad>(point_t(min.x(), min.y(), max.z()),  dx,  dy, material)); // front
        push_back(std::make_shared<Quad>(point_t(max.x(), min.y(), max.z()), -dz,  dy, material)); // right
        push_back(std::make_shared<Quad>(point_t(max.x(), min.y(), min.z()), -dx,  dy, material)); // back
        push_back(std::make_shared<Quad>(point_t(min.x(), min.y(), min.z()),  dz,  dy, material)); // left
        push_back(std::make_shared<Quad>(point_t(min.x(), max.y(), max.z()),  dx, -dz, material)); // top
        push_back(std::make_shared<Quad>(point_t(min.x(), min.y(), min.z()),  dx,  dz, material)); // bottom
    }
};

#endif