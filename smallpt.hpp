#include "include/rtw_stb_image.h"
#include "include/load_mesh.hpp"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cassert>
#include <iostream>
#include <memory>
#include <limits>
#include <vector>
#include <algorithm>


const double INF = std::numeric_limits<double>::infinity(); // for double
const int ray_tracing_depth = 5;
const double rl_init = 0.98; // the initial roussian roulette weight
const double cost_trav = 0.125; // the cost of traversal of SAH
const double cost_isect = 1.0; // the cost of ray intersecting an object in SAH
const int n_bucket = 12; // the number of buckets in SAH
const double mis_beta = 2.0; // the power of balance heuristics in MIS
const double height_coef = 6.0; // convert rgb greyscale to actual height (a proportional height mapping)

inline double degrees_to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}
inline double random_double() {
    return std::rand() / (RAND_MAX + 1.0);
}
inline double random_double(double min, double max) {
    return min + (max-min)*random_double();
}
// [min, max]
inline int random_int(int min, int max) {
    return int(random_double(min, max+1));
}
inline double linear_to_gamma(double linear_component)
{
    if (linear_component > 0)
        return std::sqrt(linear_component);

    return 0;
}
inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

struct Vec {        // Usage: time ./explicit 16 && xv image.ppm
  double x, y, z;                  // position, also color (r,g,b)
  Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
  Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
  Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
  Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
  Vec operator/(double b) const {assert(b!=0); return Vec(x/b,y/b,z/b);}
  Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
  Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
  double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; }
  double len() const {return std::sqrt(x*x+y*y+z*z);}
  Vec operator%(const Vec&b) const {return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);} // cross product
  void print() const {std::clog<<"["<<x<<", "<<y<<", "<<z<<"]";}
};

inline Vec unit_vec(const Vec& v) {
    return v / v.len();
}

inline Vec rand_unit_vec() {
    while (true) {
        auto p = Vec(random_double(-1, 1), random_double(-1, 1), random_double(-1, 1));
        auto len_sq = p.dot(p);
        if (1e-160 < len_sq && len_sq <= 1.0)
            return p / std::sqrt(len_sq);
    }
}

inline Vec rand_cosine_direction() {
    auto r1 = random_double();
    auto r2 = random_double();

    auto phi = 2*M_PI*r1;
    auto x = std::cos(phi) * std::sqrt(r2);
    auto y = std::sin(phi) * std::sqrt(r2);
    auto z = std::sqrt(1-r2);

    return Vec(x, y, z);
}

inline Vec reflect(const Vec& in, const Vec& n) {
    return in - n * ((in.dot(n)) * 2);
}

inline Vec refract(const Vec& in, const Vec& n, double r) {
    double cos_theta = std::fmin(-in.dot(n), 1.0);
    Vec out_perp = (in + n * cos_theta) * r;
    Vec out_para = n * std::sqrt(std::fabs(1.0 - out_perp.dot(out_perp)));
    return out_perp - out_para;
}

inline Vec tungsten2vec(const Tungsten::Vec3f& tv) {
    return Vec(tv.x(), tv.y(), tv.z());
}

struct Ray { 
    Vec o, d; // origin, direction
    double t; // time [0, 1]

    Ray() {}
    Ray(const Vec& o_, const Vec& d_) : o(o_), d(d_),  t(0) {} 
    Ray(const Vec& o_, const Vec& d_, double t_): o(o_), d(d_),  t(t_) {} 

    // Return: o + dist * d
    Vec at(double dist) const {
        return o + d * dist;
    }
};

struct Texture;
struct Material;
struct HitRecord {
    Vec ip, n; // intersection_point, normal (reverse direction as r_in; always normalized)
    Vec t_u, b_u; // TBN matrix. The tangent space. Note that they are normalized. They are used to convert normal map in the tangent space to the world space.
    double u, v; // coordinate at the texture space
    bool front_face; // determine whether the object emits light
    double dist; // ip = ray.orig + n * dist
    std::shared_ptr<Material> mat;

    // Set `n` such that `dot(n, r_in) >= 0`
    void set_norm(const Ray& r_in, const Vec& out_norm_u);

    void print() {
        std::clog<<"HitRecord:\n\tintersection: ";
        ip.print();
    }
};

struct Texture {
    virtual ~Texture() = default;
    virtual Vec at(double u, double v, const Vec& p) const = 0;
    virtual Vec norm(double u, double v) const = 0;
    virtual double height(double u, double v) const = 0;
};

struct ColorTexture : Texture {
    Vec color;
    ColorTexture(const Vec& c_): color(c_) {}

    // Return the color at (u, v)
    Vec at(double u, double v, const Vec& p) const override { return color; }
    Vec norm(double u, double v) const override { return Vec(); }
    double height(double u, double v) const override { return 0.0; }
};

// Copied from the textbook
struct ImageTexture : Texture {
    rtw_image img;
    rtw_image normal_map;
    rtw_image height_map;
    ImageTexture(const char* img_path): img(img_path) {}
    ImageTexture(const char* img_path, const char* normal_path): img(img_path), normal_map(normal_path) {}
    ImageTexture(const char* img_path, const char* normal_path, const char* height_path): img(img_path), normal_map(normal_path), height_map(height_path) {}
    Vec at(double u, double v, const Vec& p) const override {
        if(img.height() <=0) {return Vec(0, 1, 1);}
        u = clamp(u, 0.0, 1.0);
        v = 1 - clamp(v, 0.0, 1.0);
        auto pixel = img.pixel_data(int(u * img.width()), int(v * img.height()));
        return Vec(pixel[0] / 255.0, pixel[1] / 255.0, pixel[2] / 255.0);
    }
    // Return the non-unit normal at (u, v) in the tangent space
    Vec norm(double u, double v) const override {
        if (normal_map.is_empty()) {return Vec();}
        // Convert RGB to coordinates in the tangent space: a unit vector
        auto pixel = normal_map.pixel_data(int(u * normal_map.width()), int(v * normal_map.height()));
        return Vec((pixel[0] * 2.0 / 255.0) - 1.0, (pixel[1] * 2.0 / 255.0) - 1.0, (pixel[2] * 2.0 / 255.0) - 1.0);
    }
    // Return the height at (u, v) in the tangent space
    double height(double u, double v) const override {
        if (height_map.is_empty()) {return 0.0;}
        // Convert RGB to coordinates in the tangent space: a unit vector. White means higher.
        auto pixel = height_map.pixel_data(int(u * height_map.width()), int(v * height_map.height()));
        double avg_h = (pixel[0] + pixel[1] + pixel[2]) / (3.0 * 255.0);
        return height_coef * std::fmin(1.0, std::fmax(0.0, avg_h));
    }
};

// BSDF framework
// `sample`: sampling an outward direction for the incident ray
// `pScatter`: caculate the probability of (r_in, r_out)
struct Material {
    virtual ~Material() = default;
    virtual double pScatter(const Ray& r_in, const Ray& r_out, const HitRecord& hit_rec) const = 0;
    virtual Vec emit(const HitRecord& hit_rec) const = 0;
    virtual Vec attenuation(const HitRecord& hit_rec) const = 0;
    virtual Vec sample(const Ray& r_in, const HitRecord& hit_rec) const = 0;
    virtual double prob(const HitRecord& hit_rec, const Vec& direction) const = 0;
    virtual double sample_thresh() const = 0;
    virtual std::shared_ptr<Texture> get_tex() const = 0;
};

void HitRecord::set_norm(const Ray& r_in, const Vec& out_norm_u) {
    // Normal mapping & Height mapping
    front_face = r_in.d.dot(out_norm_u) < 0;
    n = (front_face) ? out_norm_u : out_norm_u * (-1);
    if (mat == nullptr) {return;}
    std::shared_ptr<Texture> tex = mat->get_tex();
    if (tex == nullptr) {return;}
    // Normal mapping
    Vec tangent_norm = tex->norm(u, v);
    if (std::fabs(tangent_norm.len()) > 1e-3) {
        Vec world_norm_u = unit_vec(t_u * tangent_norm.x + b_u * tangent_norm.y + out_norm_u * tangent_norm.z);
        n = (front_face) ? world_norm_u : world_norm_u * (-1);
    }
    // Height mapping
    double height = tex->height(u, v);
    ip = ip + n * height;
}

struct DiffusiveMaterial : Material {
    std::shared_ptr<Texture> tex;
    DiffusiveMaterial(const Vec& c_): tex(std::make_shared<ColorTexture>(c_)) {}
    DiffusiveMaterial(const char* img_path): tex(std::make_shared<ImageTexture>(img_path)) {}
    DiffusiveMaterial(const char* img_path, const char* normal_path): tex(std::make_shared<ImageTexture>(img_path, normal_path)) {}
    DiffusiveMaterial(const char* img_path, const char* normal_path, const char* height_path): tex(std::make_shared<ImageTexture>(img_path, normal_path, height_path)) {}
    double pScatter(const Ray& r_in, const Ray& r_out, const HitRecord& hit_rec) const override {
        double p_scatter = hit_rec.n.dot(unit_vec(r_out.d)) / M_PI;
        return p_scatter < 0 ? 0 : p_scatter;
    }
    Vec emit(const HitRecord& hit_rec) const override {
        return Vec();
    }
    Vec attenuation(const HitRecord& hit_rec) const override {
        return tex->at(hit_rec.u, hit_rec.v, hit_rec.ip);
    }   
    Vec sample(const Ray& r_in, const HitRecord& hit_rec) const override {
        double phi = 2 * M_PI * random_double();
        double r = random_double();
        double x = std::cos(phi) * std::sqrt(r);
        double y = std::sin(phi) * std::sqrt(r);
        double z = std::sqrt(1 - r);
        Vec w = unit_vec(hit_rec.n);
        Vec a = (std::fabs(w.x) > 0.9) ? Vec(0, 1, 0) : Vec(1, 0, 0);
        Vec v = unit_vec(w % a);
        Vec u = w % v;
        return u * x + v * y + w * z;
    }
    double prob(const HitRecord& hit_rec, const Vec& direction) const override {
        Vec w = unit_vec(hit_rec.n);
        return std::fmax(0, w.dot(unit_vec(direction))) / M_PI;
    }
    double sample_thresh() const override {return 0.5;}
    std::shared_ptr<Texture> get_tex() const override {return tex;}
};

struct SpecularMaterial : Material {
    std::shared_ptr<Texture> tex;
    double fuzz;
    SpecularMaterial(const Vec& c_): tex(std::make_shared<ColorTexture>(c_)), fuzz(0.5) {}
    SpecularMaterial(const char* img_path): tex(std::make_shared<ImageTexture>(img_path)) {}
    SpecularMaterial(const char* img_path, const char* normal_path): tex(std::make_shared<ImageTexture>(img_path, normal_path)) {}
    SpecularMaterial(const char* img_path, const char* normal_path, const char* height_path): tex(std::make_shared<ImageTexture>(img_path, normal_path, height_path)) {}
    SpecularMaterial(const Vec& c_, double f_): tex(std::make_shared<ColorTexture>(c_)), fuzz(std::fmax(f_, 0.01)) {}
    double pScatter(const Ray& r_in, const Ray& r_out, const HitRecord& hit_rec) const override {return 1.0;}
    Vec emit(const HitRecord& hit_rec) const override {return Vec();}
    Vec attenuation(const HitRecord& hit_rec) const override {
        return tex->at(hit_rec.u, hit_rec.v, hit_rec.ip);
    }
    Vec sample(const Ray& r_in, const HitRecord& hit_rec) const override {
        return unit_vec(reflect(r_in.d, hit_rec.n)) + rand_unit_vec() * fuzz;
    }
    double prob(const HitRecord& hit_rec, const Vec& direction) const override {return 1.0;}
    double sample_thresh() const override {return 0.0;}
    std::shared_ptr<Texture> get_tex() const override {return tex;}
};

struct TransmissiveMaterial : Material {
    double refr;
    TransmissiveMaterial(): refr(0.5) {}
    TransmissiveMaterial(double r_): refr(r_) {}
    double pScatter(const Ray& r_in, const Ray& r_out, const HitRecord& hit_rec) const override {return 1.0;}
    Vec emit(const HitRecord& hit_rec) const override {return Vec();}
    Vec attenuation(const HitRecord& hit_rec) const override {return Vec(1.0, 1.0, 1.0);}   
    Vec sample(const Ray& r_in, const HitRecord& hit_rec) const override {
        double refr_factor = hit_rec.front_face ? (1.0 / refr) : refr;
        Vec in_dir_u = unit_vec(r_in.d);
        double cos_theta = std::fmin(-in_dir_u.dot(hit_rec.n), 1.0);
        double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
        if (refr_factor * sin_theta > 1.0 || schlick_approx(cos_theta, refr_factor) > random_double()) {
            return reflect(in_dir_u, hit_rec.n);
        }
        else {
            return refract(in_dir_u, hit_rec.n, refr_factor);
        }
    }
    double prob(const HitRecord& hit_rec, const Vec& direction) const override {return 1.0;}
    double sample_thresh() const override {return 0.0;}
    std::shared_ptr<Texture> get_tex() const override {return nullptr;}
    double schlick_approx(double cos_theta, double refr_factor) const {
        double r0 = (1 - refr_factor) / (1 + refr_factor);
        r0 *= r0;
        return r0 + (1 - r0) * std::pow((1 - cos_theta), 5);
    }
};

struct AreaLight : Material {
    std::shared_ptr<Texture> tex;
    AreaLight(const Vec& c_): tex(std::make_shared<ColorTexture>(c_)) {}
    AreaLight(const char* img_path): tex(std::make_shared<ImageTexture>(img_path)) {}
    double pScatter(const Ray& r_in, const Ray& r_out, const HitRecord& hit_rec) const override {assert(false); return 0.0;}
    Vec emit(const HitRecord& hit_rec) const override {
        if (!hit_rec.front_face){return Vec();}
        return tex->at(hit_rec.u, hit_rec.v, hit_rec.ip);
    }
    Vec attenuation(const HitRecord& hit_rec) const override {assert(false);return Vec(1.0, 1.0, 1.0);}  
    Vec sample(const Ray& r_in, const HitRecord& hit_rec) const override {return Vec();} 
    double prob(const HitRecord& hit_rec, const Vec& direction) const override {return 0.0;}
    double sample_thresh() const override {return 0.0;}
    std::shared_ptr<Texture> get_tex() const override {return tex;}
};

struct BBox {
    double x0, x1, y0, y1, z0, z1;

    BBox(double x0_=INF, double x1_=-INF, double y0_=INF, double y1_=-INF, double z0_=INF, double z1_=-INF) {
        set(x0_, x1_, y0_, y1_, z0_, z1_);
    }

    BBox(const Vec& a, const Vec& b) {
        double x0_ = a.x <= b.x ? a.x : b.x;
        double x1_ = a.x > b.x ? a.x : b.x;
        double y0_ = a.y <= b.y ? a.y : b.y;
        double y1_ = a.y > b.y ? a.y : b.y;
        double z0_ = a.z <= b.z ? a.z : b.z;
        double z1_ = a.z > b.z ? a.z : b.z;
        set(x0_, x1_, y0_, y1_, z0_, z1_);
    }

    BBox(const BBox& a, const BBox& b) {
        double x0_ = a.x0 <= b.x0 ? a.x0 : b.x0;
        double x1_ = a.x1 > b.x1 ? a.x1 : b.x1;
        double y0_ = a.y0 <= b.y0 ? a.y0 : b.y0;
        double y1_ = a.y1 > b.y1 ? a.y1 : b.y1;
        double z0_ = a.z0 <= b.z0 ? a.z0 : b.z0;
        double z1_ = a.z1 > b.z1 ? a.z1 : b.z1;
        set(x0_, x1_, y0_, y1_, z0_, z1_);
    }
    
    // Return: whether `r` intersect this bounding box
    bool intersect(const Ray& r, double min, double max) const {
        // x
        double x_range_0 = (x0 - r.o.x) / r.d.x;
        double x_range_1 = (x1 - r.o.x) / r.d.x;
        double x_range_min = x_range_0 <= x_range_1 ? x_range_0 : x_range_1;
        double x_range_max = x_range_0 > x_range_1 ? x_range_0 : x_range_1;
        // std::clog<<"\t"<<x_range_0<<" "<<x_range_1<<"\n";/////////////////

        // y
        double y_range_0 = (y0 - r.o.y) / r.d.y;
        double y_range_1 = (y1 - r.o.y) / r.d.y;
        double y_range_min = y_range_0 <= y_range_1 ? y_range_0 : y_range_1;
        double y_range_max = y_range_0 > y_range_1 ? y_range_0 : y_range_1;
        // std::clog<<"\t"<<y_range_min<<" "<<y_range_max<<"\n";/////////////////

        // z
        double z_range_0 = (z0 - r.o.z) / r.d.z;
        double z_range_1 = (z1 - r.o.z) / r.d.z;
        double z_range_min = z_range_0 <= z_range_1 ? z_range_0 : z_range_1;
        double z_range_max = z_range_0 > z_range_1 ? z_range_0 : z_range_1;
        // std::clog<<"\t"<<z_range_min<<" "<<z_range_max<<"\n";/////////////////

        double range_min_max = (x_range_min > y_range_min && x_range_min > z_range_min) ? x_range_min : (y_range_min > z_range_min ? y_range_min : z_range_min);
        double range_max_min = (x_range_max < y_range_max && x_range_max < z_range_max) ? x_range_max : (y_range_max < z_range_max ? y_range_max : z_range_max);
        // std::clog<<range_min_max<<" "<<range_max_min<<std::endl;////////////
        range_min_max = range_min_max > min ? range_min_max : min;
        range_max_min = range_max_min < max ? range_max_min : max;
        return range_min_max < range_max_min;
    }

    int longest_axis() const {
        double range_x = x1 - x0;
        double range_y = y1 - y0;
        double range_z = z1 - z0;
        return (range_x > range_y && range_x > range_z) ? 0 : (range_y > range_z ? 1 : 2);
    }

    // For SAH
    double surface_area() const {
        double area = std::fmax(x1 - x0, 0) * std::fmax(y1 - y0, 0) + std::fmax(x1 - x0, 0) * std::fmax(z1 - z0, 0) + std::fmax(y1 - y0, 0) * std::fmax(z1 - z0, 0);
        assert(area >= 0);
        return 2.0 * area;
    }
    // For SAH
    double centroid(int axis) const {
        switch (axis)
        {
        case 0:
            return (x1 + x0) / 2.0;
        case 1:
            return (y1 + y0) / 2.0;
        case 2:
            return (z1 + z0) / 2.0;
        default:
            std::cerr << "Unknown axis " << axis << std::endl;
            exit(1);
        }
    }

    void print() const {
        std::clog<<"["<<x0<<", "<<x1<<"] x ";
        std::clog<<"["<<y0<<", "<<y1<<"] x ";
        std::clog<<"["<<z0<<", "<<z1<<"]";
    }

    void set(double x0_, double x1_, double y0_, double y1_, double z0_, double z1_) {
        x0=x0_; x1=x1_; y0=y0_; y1=y1_; z0=z0_; z1=z1_;
        double delta = 0.01;
        if (x1 - x0 < delta) {x0 -= delta/2; x1 += delta/2;}
        if (y1 - y0 < delta) {y0 -= delta/2; y1 += delta/2;}
        if (z1 - z0 < delta) {z0 -= delta/2; z1 += delta/2;}
    }
};

struct Object {
    virtual ~Object() = default;
    virtual bool intersect(const Ray &r, HitRecord& hit_rec, double min, double max) const = 0;
    virtual double prob(const Vec& src, const Vec& direction) const = 0; 
    virtual Vec sample(const Vec& src) const = 0;
    virtual BBox get_bbox() const = 0;
};

struct Sphere : Object {
    double rad;       // radius
    Ray o;      // center. Support motion.
    std::shared_ptr<Material> mat;      // reflection type (DIFFuse, SPECular, REFRactive)
    BBox bbox;
    double is_hdr; // whether this is the HDR World sphere

    // Static sphere
    Sphere(const Vec& static_center, double radius, std::shared_ptr<Material> material) : o(Ray(static_center, Vec())), rad(std::fmax(radius, 0)), mat(material), is_hdr(false) {
        Vec radius_v = Vec(rad, rad, rad);
        bbox = BBox(static_center - radius_v, static_center + radius_v);
    }

    // Moving sphere
    Sphere(const Vec& from, const Vec& to, double radius, std::shared_ptr<Material> material) : o(Ray(from, to - from)), rad(std::fmax(radius, 0)), mat(material), is_hdr(false) {
        Vec radius_v = Vec(rad, rad, rad);
        bbox = BBox(BBox(from - radius_v, from + radius_v), BBox(to - radius_v, to + radius_v));
    }

    // HDR Static sphere
    Sphere(const Vec& static_center, double radius, std::shared_ptr<Material> material, bool hdr_) : o(Ray(static_center, Vec())), rad(std::fmax(radius, 0)), mat(material), is_hdr(hdr_) {
        Vec radius_v = Vec(rad, rad, rad);
        bbox = BBox(static_center - radius_v, static_center + radius_v);
    }

    // Return whether intersection exists
    // Store `hit_rec`
    bool intersect(const Ray &r, HitRecord& hit_rec, double min, double max) const override { // returns distance, 0 if nohit
        Vec cur_o = o.at(r.t);
        Vec ray2sphere = cur_o - r.o;
        double b = ray2sphere.dot(r.d);
        double a = r.d.dot(r.d);
        double c = ray2sphere.dot(ray2sphere) - rad * rad;
        double det = b * b - a * c;
        // double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
        // if (det<0) return 0; else det=sqrt(det);
        // return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
        if (det < 0) {return false;}
        double sqrt_det = std::sqrt(det);
        double root = (b - sqrt_det) / a;
        if (root <= min || root > max) {root = (b + sqrt_det) / a;}
        if (root <= min || root > max) {return false;}
        Vec ip = r.at(root);
        Vec out_norm_u = (ip - cur_o) / rad;
        double theta = std::acos(-out_norm_u.y);
        double phi = std::atan2(-out_norm_u.z, out_norm_u.x) + M_PI;
        double cos_theta = std::cos(theta);
        double sin_theta = std::sin(theta);
        double cos_phi = std::cos(phi);
        double sin_phi = std::sin(phi);
        hit_rec.u = phi / (2 * M_PI);
        hit_rec.v = theta / M_PI;
        hit_rec.t_u = Vec(sin_phi, 0, cos_phi);
        hit_rec.b_u = Vec(-sin_theta * cos_phi, cos_theta, sin_theta * sin_phi);
        hit_rec.dist = root;
        hit_rec.ip = ip;
        hit_rec.mat = mat;
        if (!is_hdr) {hit_rec.set_norm(r, out_norm_u);}
        else {hit_rec.set_norm(r, out_norm_u * (-1));}
        return true;
    }

    // Importance sampling; calc the pdf for lights of various shapes
    // Return: the prob of getting `direction`
    // Only for static sphere
    double prob(const Vec& src, const Vec& direction) const override {
        // assert(o.d.len() == 0);
        HitRecord hit_rec;
        if (!this->intersect(Ray(src, direction), hit_rec, 1e-3, INF)) {
            return 0;
        }
        if (is_hdr) {
            return 1 / (4 * M_PI);
        }
        // not HDR sphere
        Vec src2sphere = o.at(0) - src;
        double cos_theta = std::sqrt(1 - rad * rad / (src2sphere.dot(src2sphere)));
        return 1 / (2 * M_PI * (1 - cos_theta));
    }

    // Importance sampling; sample a point on the shape
    // Return: a sampled `direction`
    // Only for static sphere
    Vec sample(const Vec& src) const override {
        if (is_hdr) {
            // For HDR, the point `src` is gauaranteed to be inside this sphere. 
            // So, we uniformly sample a point on the sphere.
            // A more sophisticated HDR light uses mipmap for sampling. Currently I haven't implemented this yet.
            Vec rand_sphere_point = o.at(0) + rand_unit_vec() * rad;
            return rand_sphere_point - src;
        }
        // not HDR sphere
        Vec src2sphere = o.at(0) - src;
        double dist_sq = src2sphere.dot(src2sphere);
        double z = 1 + random_double() * (std::sqrt(1 - rad * rad / dist_sq) - 1);
        double phi = 2 * M_PI * random_double();
        double x = std::cos(phi) * std::sqrt(1 - z * z);
        double y = std::sin(phi) * std::sqrt(1 - z * z);
        Vec w = unit_vec(src2sphere);
        Vec a = (std::fabs(w.x) > 0.9) ? Vec(0, 1, 0) : Vec(1, 0, 0);
        Vec v = unit_vec(w % a);
        Vec u = w % v;
        return u * x + v * y + w * z;
    }

    BBox get_bbox() const override {return bbox;}

};

struct Quad : Object {
    Vec o, u, v;
    Vec u_u, v_u;
    Vec uv_cross; // `n_u` is the unit-len normal following right-hand rule
    double dot_o_n_u, uv_cross_len_sq; // <n_u, o>
    Vec n_u; // face normal of unit len; follow right-hand rule
    std::shared_ptr<Material> mat;
    BBox bbox;

    Quad(const Vec& o_, const Vec& u_, const Vec& v_, std::shared_ptr<Material> mat_) : o(o_), u(u_), v(v_), mat(mat_) {
        uv_cross = u % v;
        uv_cross_len_sq = uv_cross.dot(uv_cross);
        n_u = unit_vec(uv_cross);
        dot_o_n_u = n_u.dot(o);
        bbox = BBox(BBox(o, o+u+v), BBox(o+u, o+v));
        u_u = unit_vec(u);
        v_u = unit_vec(v);
    }

    // Return whether intersection exists
    // Store `hit_rec`
    bool intersect(const Ray &r, HitRecord& hit_rec, double min, double max) const override {
        double dot_r_n_u = r.d.dot(n_u);
        if(std::fabs(dot_r_n_u) < 1e-8) {return false;} // parallel
        double dist = (dot_o_n_u - n_u.dot(r.o)) / dot_r_n_u;
        if(dist <= min || dist >= max) {return false;}
        Vec ip = r.at(dist);
        Vec ip_offset = ip - o;
        double alpha = (uv_cross.dot(ip_offset % v)) / uv_cross_len_sq;
        double beta = (uv_cross.dot(u % ip_offset)) / uv_cross_len_sq;
        if(alpha>1 || alpha<0 || beta>1 || beta<0) {return false;}
        hit_rec.u = alpha;
        hit_rec.v = beta;
        hit_rec.t_u = u_u;
        hit_rec.b_u = v_u;
        hit_rec.dist = dist;
        hit_rec.ip = ip;
        hit_rec.mat = mat;
        hit_rec.set_norm(r, n_u);
        return true;
    }

    // Importance sampling; generate the pdf for lights of various shapes
    // Return: the probability of getting `direction`
    double prob(const Vec& src, const Vec& direction) const override {
        HitRecord hit_rec;
        if (!this->intersect(Ray(src, direction), hit_rec, 1e-3, INF)) {
            return 0;
        }
        double area = std::sqrt(uv_cross_len_sq);
        double direction_len = std::sqrt(direction.dot(direction));
        return (hit_rec.dist * hit_rec.dist * direction_len * direction_len * direction_len) / (area * std::fabs(direction.dot(hit_rec.n)));
    }

    // Importance sampling; sample a point on the shape
    // Return: a sampled `direction`
    Vec sample(const Vec& src) const override {
        return o + u * random_double() + v * random_double() - src;
    }
    
    BBox get_bbox() const override {return bbox;}

};

struct Transform: Object {
    Vec trans;
    std::shared_ptr<Object> obj;
    BBox bbox;
    Transform(const Vec& trans_, std::shared_ptr<Object> obj_): trans(trans_), obj(obj_) {
        bbox = obj->get_bbox();
        bbox.x0 += trans.x;
        bbox.y0 += trans.y;
        bbox.z0 += trans.z;
        bbox.x1 += trans.x;
        bbox.y1 += trans.y;
        bbox.z1 += trans.z;
    }
    bool intersect(const Ray &r, HitRecord& hit_rec, double min, double max) const override {
        Ray offset_ray(r.o - trans, r.d, r.t);
        if (!obj->intersect(offset_ray, hit_rec, min, max)) {return false;}
        hit_rec.ip = hit_rec.ip + trans;
        return true;
    }
    double prob(const Vec& src, const Vec& direction) const override {
        Vec offset_src = src - trans;
        return obj->prob(offset_src, direction);
    }; 
    Vec sample(const Vec& src) const override {
        Vec offset_src = src - trans;
        return obj->sample(offset_src);
    };
    BBox get_bbox() const override {return bbox;};
};

struct Triangle : Object {
    std::shared_ptr<Tungsten::Vertex> o;
    std::shared_ptr<Tungsten::Vertex> v_u;
    std::shared_ptr<Tungsten::Vertex> v_v;
    Vec uv_cross;
    double dot_o_n_u, uv_cross_len_sq;
    Vec n_u;
    std::shared_ptr<Material> mat;
    BBox bbox;
    bool is_valid;

    Triangle(const std::vector<Tungsten::Vertex>& vtab, const Tungsten::TriangleI& tri, std::shared_ptr<Material> mat_): mat(mat_) {
        o = std::make_shared<Tungsten::Vertex>(vtab[tri.v0]);
        Vec o_norm = tungsten2vec(o->normal());
        Vec o_pos = tungsten2vec(o->pos());
        Vec v1 = tungsten2vec(vtab[tri.v1].pos()) - o_pos;
        Vec v2 = tungsten2vec(vtab[tri.v2].pos()) - o_pos;
        uv_cross = v1 % v2;
        if (uv_cross.dot(o_norm) < 0) {
            uv_cross = uv_cross * (-1);
            v_v = std::make_shared<Tungsten::Vertex>(vtab[tri.v1]);
            v_u = std::make_shared<Tungsten::Vertex>(vtab[tri.v2]);
        } else {
            v_u = std::make_shared<Tungsten::Vertex>(vtab[tri.v1]);
            v_v = std::make_shared<Tungsten::Vertex>(vtab[tri.v2]);
        }
        uv_cross_len_sq = uv_cross.dot(uv_cross);
        if (uv_cross_len_sq == 0) {is_valid = false; return;}
        else {is_valid = true;}
        n_u = unit_vec(uv_cross);
        dot_o_n_u = n_u.dot(o_pos);
        bbox = BBox(BBox(o_pos, o_pos + v1 + v2), BBox(o_pos + v1, o_pos + v2));
        
    }

    bool intersect(const Ray &r, HitRecord& hit_rec, double min, double max) const override {
        Vec o_pos = tungsten2vec(o->pos());
        Vec u_dir = tungsten2vec(v_u->pos()) - o_pos;
        Vec v_dir = tungsten2vec(v_v->pos()) - o_pos;

        double dot_r_n_u = r.d.dot(n_u);
        if(std::fabs(dot_r_n_u) < 1e-8) {return false;} // parallel
        double dist = (dot_o_n_u - n_u.dot(r.o)) / dot_r_n_u;
        if(dist <= min || dist >= max) {return false;}
        Vec ip = r.at(dist);
        Vec ip_offset = ip - o_pos;
        double alpha = (uv_cross.dot(ip_offset % v_dir)) / uv_cross_len_sq;
        double beta = (uv_cross.dot(u_dir % ip_offset)) / uv_cross_len_sq;
        if(alpha>1 || alpha<0 || beta>1 || beta<0 || (alpha + beta)>1) {return false;}
        hit_rec.u = alpha;
        hit_rec.v = beta;
        hit_rec.t_u = unit_vec(u_dir);
        hit_rec.b_u = unit_vec(v_dir);
        hit_rec.dist = dist;
        hit_rec.ip = ip;
        hit_rec.mat = mat;
        hit_rec.set_norm(r, n_u);
        return true;
    };

    double prob(const Vec& src, const Vec& direction) const override {
        HitRecord hit_rec;
        if (!this->intersect(Ray(src, direction), hit_rec, 1e-3, INF)) {
            return 0;
        }
        double area = std::sqrt(uv_cross_len_sq) / 2.0;
        double direction_len = std::sqrt(direction.dot(direction));
        return (hit_rec.dist * hit_rec.dist * direction_len * direction_len * direction_len) / (area * std::fabs(direction.dot(hit_rec.n)));
    }; 

    Vec sample(const Vec& src) const override {
        Vec o_pos = tungsten2vec(o->pos());
        Vec u_dir = tungsten2vec(v_u->pos()) - o_pos;
        Vec v_dir = tungsten2vec(v_v->pos()) - o_pos;

        return o_pos + (u_dir * random_double() + v_dir * random_double()) * 0.5 - src;
    };

    BBox get_bbox() const override {return bbox;}
};

void create_mesh(std::vector<std::shared_ptr<Object>>& objects, const std::string& relative_mesh_wo3, std::vector<Tungsten::Vertex>& vtab, std::vector<Tungsten::TriangleI>& ftab, const std::shared_ptr<Material> material, const Vec& trans) {
    bool succeed_loading_mesh = TungstenloadWo3(relative_mesh_wo3, vtab, ftab);
    if (!succeed_loading_mesh) {
        std::clog <<"Fail to load the mesh from "<<relative_mesh_wo3<<"\n";
        return;
    }
    double trans_len = trans.len();
    for(auto& tri : ftab) {
        auto my_tri = std::make_shared<Triangle>(vtab, tri, material);
        if (my_tri->is_valid == true && trans_len == 0.0) {
            objects.push_back(my_tri);
        } 
        else if (my_tri->is_valid == true) {
            objects.push_back(std::make_shared<Transform>(trans, my_tri));
        }
    }

}

void create_box(std::vector<std::shared_ptr<Object>>& objects, const Vec& a, double x, double y, double z, const std::shared_ptr<Material> material) {
    Vec b = a + Vec(std::fabs(x), std::fabs(y), std::fabs(z));
    ////////////////
    // a.print();
    // std::clog<<"\n";
    // b.print();
    //////////////////
    Vec min = Vec(std::fmin(a.x,b.x), std::fmin(a.y,b.y), std::fmin(a.z,b.z));
    Vec max = Vec(std::fmax(a.x,b.x), std::fmax(a.y,b.y), std::fmax(a.z,b.z));

    Vec dx = Vec(max.x - min.x, 0, 0);
    Vec dy = Vec(0, max.y - min.y, 0);
    Vec dz = Vec(0, 0, max.z - min.z);

    objects.push_back(std::make_shared<Quad>(Vec(min.x, min.y, max.z),  dx,  dy, material)); // front
    objects.push_back(std::make_shared<Quad>(Vec(max.x, min.y, max.z), dz * (-1),  dy, material)); // right
    objects.push_back(std::make_shared<Quad>(Vec(max.x, min.y, min.z), dx * (-1),  dy, material)); // back
    objects.push_back(std::make_shared<Quad>(Vec(min.x, min.y, min.z),  dz,  dy, material)); // left
    objects.push_back(std::make_shared<Quad>(Vec(min.x, max.y, max.z),  dx, dz * (-1), material)); // top
    objects.push_back(std::make_shared<Quad>(Vec(min.x, min.y, min.z),  dx,  dz, material)); // bottom  
}

struct BucketInfo {
    int n_objects;
    BBox bbox;
    int obj_idx;

    BucketInfo(): n_objects(0), bbox(BBox()), obj_idx(-1) {}
};

// BVH with SAH (surface area heuristic)
struct BVHNode {
    int start, end; // idx range of covered objects
    std::shared_ptr<BVHNode> left;
    std::shared_ptr<BVHNode> right;
    BBox bbox;
    double is_bvh;

    BVHNode(std::vector<std::shared_ptr<Object>>& objects, int start_idx, int end_idx): start(start_idx), end(end_idx), is_bvh(end_idx - start_idx > 1) {
        assert(end_idx > start_idx);
        // Set bbox
        bbox = BBox();
        for(int idx = start_idx; idx < end_idx; idx++) {
            bbox = BBox(bbox, objects[idx]->get_bbox());
        }
        
        // Recursively build children
        // base
        if (!is_bvh) {
            assert(end_idx - start_idx == 1);
            return;
        }
        /////////////////////////////
        // std::clog<<"Create BBox for bvh node:\n\t";/////////////////////////////
        // bbox.print();//////////////////////////////////////////////////////
        // std::clog<<"\n";
        //////////////////////////////////
        assert(end_idx - start_idx >= 2);
        int max_axis = bbox.longest_axis();
        switch (max_axis)
        {
        case 0: //x
            std::sort(std::begin(objects) + start_idx, std::begin(objects) + end_idx, box_x_compare);
            break;
        case 1: // y
            std::sort(std::begin(objects) + start_idx, std::begin(objects) + end_idx, box_y_compare);
            break;
        case 2: // z
            std::sort(std::begin(objects) + start_idx, std::begin(objects) + end_idx, box_z_compare);
            break;
        default:
            std::cerr << "Unknown axis " << max_axis << std::endl;
            exit(1);
        }
        // SAH: how to choose `mid`
        int mid = start_idx + (end_idx - start_idx) / 2; // default: end_idx - start_idx <= 4
        if (end_idx - start_idx > 4) {
            // Create buckets
            double left_margin = objects[start_idx]->get_bbox().centroid(max_axis);
            double right_margin = objects[end_idx - 1]->get_bbox().centroid(max_axis);
            double bucket_size = (right_margin - left_margin) / n_bucket;
            if (bucket_size <= 0) {
                // all aligned
                mid = start_idx + (end_idx - start_idx) / 2;
            }
            else {
                assert(bucket_size > 0);
                //Initialize the buckets: O(n)
                std::vector<BucketInfo> buckets(n_bucket);
                for (int idx = start_idx; idx < end_idx; idx++) {
                    double centroid = objects[idx]->get_bbox().centroid(max_axis);
                    // Predict the bucket idx
                    int bucket_idx = int((centroid - left_margin) / bucket_size);
                    assert(bucket_idx >= 0);
                    if (bucket_idx == n_bucket) {bucket_idx = n_bucket - 1;}
                    assert(bucket_idx < n_bucket);
                    // Update the bucket
                    buckets[bucket_idx].n_objects += 1;
                    buckets[bucket_idx].bbox = BBox(buckets[bucket_idx].bbox, objects[idx]->get_bbox());
                    buckets[bucket_idx].obj_idx = idx;
                }
                // For debugging
                for (int bid = 0; bid < n_bucket - 1; bid++) {
                    assert(buckets[bid + 1].obj_idx == -1 || buckets[bid].obj_idx == -1 || buckets[bid].obj_idx < buckets[bid + 1].obj_idx);
                }
                // Find the best `mid`: O(n)
                double tot_area = bbox.surface_area();
                assert(tot_area > 0);
                double min_cost = INF;
                double opt_mid = -1; // #left
                for (int cur_mid = 1; cur_mid < n_bucket; cur_mid++) {
                    BBox left_bbox = BBox();
                    int left_n_objects = 0;
                    int left_mid = -1;
                    BBox right_bbox = BBox();
                    int right_n_objects = 0;
                    for (int lidx = 0; lidx < cur_mid; lidx++) {
                        left_bbox = BBox(left_bbox, buckets[lidx].bbox);
                        left_n_objects += buckets[lidx].n_objects;
                        if (buckets[lidx].obj_idx != -1) {
                            assert(buckets[lidx].obj_idx > left_mid);
                            left_mid = buckets[lidx].obj_idx;
                        }
                    }
                    if(left_mid == -1) {assert(left_n_objects==0); continue;}
                    for (int ridx = cur_mid; ridx < n_bucket; ridx++) {
                        right_bbox = BBox(right_bbox, buckets[ridx].bbox);
                        right_n_objects += buckets[ridx].n_objects;
                    }
                    if(right_n_objects == 0) {continue;}
                    assert(left_n_objects + right_n_objects == end_idx - start_idx);
                    double cur_cost = cost_trav + cost_isect * (left_n_objects * left_bbox.surface_area() + right_n_objects * right_bbox.surface_area()) / tot_area;
                    if (cur_cost < min_cost) {
                        min_cost = cur_cost;
                        opt_mid = left_mid;
                    }
                }
                assert(min_cost < INF);
                assert(opt_mid >= start_idx and opt_mid <= end_idx - 2);
                mid = opt_mid + 1;
            }
        }
        left = std::make_shared<BVHNode>(objects, start_idx, mid);
        right = std::make_shared<BVHNode>(objects, mid, end_idx);
    }
    // Sort according to centroid, due to SAH
    static bool box_x_compare (const std::shared_ptr<Object> a, const std::shared_ptr<Object> b) {return a->get_bbox().centroid(0) < b->get_bbox().centroid(0);}
    static bool box_y_compare (const std::shared_ptr<Object> a, const std::shared_ptr<Object> b) {return a->get_bbox().centroid(1) < b->get_bbox().centroid(1);}
    static bool box_z_compare (const std::shared_ptr<Object> a, const std::shared_ptr<Object> b) {return a->get_bbox().centroid(2) < b->get_bbox().centroid(2);}
};

// Multiple importance sampling
// Return the MIS weight: $\frac{ N (n_i p_i) ^ {b - 1} }{ \sum_k n_kp_k }$
// Store the sampled ray to `r_out`
double sample(
    const Ray& r_in, 
    const HitRecord& hit_rec, 
    const std::vector<std::shared_ptr<Object>>& point_lights, 
    Ray& r_out,
    int obj_idx,
    int n_samples,
    const std::vector<int>& mis_samples
) {
    // Sample a scattering direction
    r_out.o = hit_rec.ip;
    r_out.t = r_in.t;
    double thresh = hit_rec.mat->sample_thresh();
    if (point_lights.size() == 0) {
        thresh = 0.0;
    }
    // if (random_double() < thresh) {
    //     // sample using lights based on their shapes: randomly choose a light
    //     int pl_idx = random_int(0, point_lights.size() - 1);
    //     r_out.d = point_lights[pl_idx]->sample(hit_rec.ip); // src = intersection
    // }
    // else {
    //     // sample using the material of the hit object
    //     r_out.d = hit_rec.mat->sample(r_in, hit_rec);
    // }
    // // Calculate the sampling probability
    // double pl_prob = 0.0;
    // if (thresh > 0) {
    //     for (const auto& pl : point_lights) {
    //         pl_prob += pl->prob(hit_rec.ip, r_out.d); // origin = intersection, direction = out_dir
    //     }
    //     pl_prob /= point_lights.size();
    // }
    // double mat_prob = 0.0;
    // if (1 - thresh > 0) {
    //     mat_prob = hit_rec.mat->prob(hit_rec, r_out.d);
    // }
    // return thresh * pl_prob + (1 - thresh) * mat_prob;

    // Calculate the MIS weight
    double mis_prob;
    if (obj_idx > 0 && thresh > 0.0) {
        // sample a light
        r_out.d = point_lights[obj_idx - 1]->sample(hit_rec.ip); 
        mis_prob = point_lights[obj_idx - 1]->prob(hit_rec.ip, r_out.d);
    } else {
        // sample the object
        r_out.d = hit_rec.mat->sample(r_in, hit_rec);
        mis_prob = hit_rec.mat->prob(hit_rec, r_out.d);
        if (thresh == 0.0) {
            return mis_prob == 0.0 ? 0.0 : 1 / mis_prob;
        }
        assert(obj_idx == 0);
    }
    if (mis_prob == 0.0) {return 0.0;} // emit
    double mw_n = double(n_samples) * std::pow(mis_samples[obj_idx] * mis_prob, mis_beta - 1);
    double mw_d = std::pow(mis_samples[0] * hit_rec.mat->prob(hit_rec, r_out.d), mis_beta); // the object
    for (int mw_d_idx=1; mw_d_idx < mis_samples.size(); mw_d_idx++) { // the lights
        mw_d += std::pow(mis_samples[mw_d_idx] * point_lights[mw_d_idx - 1]->prob(hit_rec.ip, r_out.d), mis_beta);
    } 
    assert(mw_n > 0 && mw_d > 0);
    return mw_n / mw_d;
}

// BVH optimization
// Return the idx of the intersected object
// If not exists, return -1
int intersect(const std::vector<std::shared_ptr<Object>>& objects, const std::shared_ptr<BVHNode> &bvh_root, const Ray &r, HitRecord& hit_rec, double min, double max) {
    if (!bvh_root->is_bvh) {
        bool is_intersect = objects[bvh_root->start]->intersect(r, hit_rec, min, max);
        return is_intersect ? bvh_root -> start : -1;
    }
    if (!bvh_root->bbox.intersect(r, min, max)) {return -1;}
    assert(bvh_root->end - bvh_root->start > 1);
    int l_intersect = intersect(objects, bvh_root->left, r, hit_rec, min, max);
    double new_max = l_intersect != -1 ? hit_rec.dist : max;
    int r_intersect = intersect(objects, bvh_root->right, r, hit_rec, min, new_max);
    return r_intersect != -1 ? r_intersect : l_intersect;
}

// Roussian Roulette
Vec radiance(
    const std::vector<std::shared_ptr<Object>>& point_lights,
    const std::vector<std::shared_ptr<Object>>& objects, 
    const std::shared_ptr<BVHNode> &bvh_root, 
    const Ray &r_in,
    int depth,
    Vec background, 
    double rl_weight,
    int mis_obj_idx,
    int n_samples,
    const std::vector<int>& mis_samples
) {
    // Roussian roulette
    if (depth <= 0 && random_double() > rl_weight) {return Vec();} // at least `ray_tracing_depth` bounces!// if (depth <= 0) {return Vec();}
    HitRecord hit_rec;
    int obj_idx = intersect(objects, bvh_root, r_in, hit_rec, 1e-3, INF);
    if (obj_idx == -1) {return background;}
    // Sample the scattered ray
    Ray r_out;
    double mis_weight = sample(r_in, hit_rec, point_lights, r_out, mis_obj_idx, n_samples, mis_samples);
    if (mis_weight == 0) {return hit_rec.mat->emit(hit_rec);}
    // Calculate pScatter, attenuation, emit
    double scattered_prob = hit_rec.mat->pScatter(r_in, r_out, hit_rec);
    Vec attenuation = hit_rec.mat->attenuation(hit_rec);
    Vec emit = hit_rec.mat->emit(hit_rec);
    // Reduce rl_weigth by attenutation
    double new_rl_weight = attenuation.x > attenuation.y && attenuation.x > attenuation.z ? attenuation.x : attenuation.y > attenuation.z ? attenuation.y : attenuation.z;
    new_rl_weight = std::fmax(new_rl_weight, 0.95);
    new_rl_weight *= rl_weight;
    new_rl_weight = std::fmin(new_rl_weight, rl_init);
    // Recursion
    Vec post_color = radiance(point_lights, objects, bvh_root, r_out, depth - 1, background, new_rl_weight, mis_obj_idx, n_samples, mis_samples);
    // Calculate the final color
    // return (emit + attenuation.mult(post_color) * scattered_prob / sampled_prob) / rl_weight;
    return (emit + attenuation.mult(post_color) * scattered_prob * mis_weight) / rl_weight;
}

// Write a pixel color to a line of a file: "r g b\n"
Vec get_color(const Vec& pixel_color) {
    double r = pixel_color.x;
    double g = pixel_color.y;
    double b = pixel_color.z;
    // Replace NaN components with zero.
    // Copied from the textbook. This is to remove black acnes.
    if (r != r) r = 0.0;
    if (g != g) g = 0.0;
    if (b != b) b = 0.0;
    r = linear_to_gamma(r);
    g = linear_to_gamma(g);
    b = linear_to_gamma(b);
    double min = 0;
    double max = 0.999;
    int rbyte = int(256 * clamp(r, min, max));
    int gbyte = int(256 * clamp(g, min, max));
    int bbyte = int(256 * clamp(b, min, max));
    return Vec(rbyte, gbyte, bbyte);
}

// Set the number of samples for each PDF
void MIS_init(std::vector<int>& mis_samples, int n_samples) {
    if (mis_samples.size() == 1) {mis_samples[0] = n_samples; return;} // there is no light
    mis_samples[0] = n_samples / 1.05;
    // mis_samples[0] = 0;
    int n_remain_samples = n_samples - mis_samples[0];
    int quotient = n_remain_samples / (mis_samples.size() - 1);
    int remainder = n_remain_samples % (mis_samples.size() - 1);
    for(int i=1; i < mis_samples.size(); i++) {
        mis_samples[i] = quotient;
        if (i <= remainder) {
            mis_samples[i] += 1;
        }
    }
    // DEBUG
    int tot_samples = 0;
    for (auto& ms : mis_samples) {
        tot_samples += ms;
    }
    assert(tot_samples == n_samples);

}

void render(
    FILE* out,
    std::vector<std::shared_ptr<Object>>& point_lights, 
    std::vector<std::shared_ptr<Object>>& objects,
    int image_height,
    int image_width,
    int sample_len,
    Vec lookfrom, // camera position
    Vec lookat, // camera target
    Vec vup, // up direction for camera
    double vfov, // camera vfov
    double focus, // focus of the camera
    Vec background_color
) {
    fprintf(out, "P3\n%d %d\n%d\n", image_width, image_height, 255);
    // Initialize the camera's view
    double view_height = std::tan(0.5 * degrees_to_radians(vfov)) * 2 * focus;
    double view_width = view_height * (double(image_width)/image_height);
    Vec w = unit_vec(lookfrom - lookat);
    Vec u = unit_vec(vup % w);
    Vec v = w % u;
    Vec view_u  = u * view_width;
    Vec view_v = v * (-view_height);
    Vec pixel_u = view_u / image_width;
    Vec pixel_v = view_v / image_height;
    Vec pixel_orig = lookfrom - w * focus - view_u / 2 - view_v / 2 + (pixel_u + pixel_v) / 2;
    // Create BVH tree for optimization
    std::shared_ptr<BVHNode> bvh_root = std::make_shared<BVHNode>(objects, 0, objects.size());
    // Raytracing
    std::vector<Vec> image(image_height * image_width);
    // Multiple importance sampling
    std::vector<int> mis_samples(1 + point_lights.size());
    MIS_init(mis_samples, sample_len * sample_len);
    std::clog<<mis_samples[0]<<"\n";////////////////////////////////////
    for(int j = 0; j < image_height; j++) {
        std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
        for(int i = 0; i < image_width; i++) {
            Vec pixel_color;
            // Anti-aliasing by sampling rays for each pixel
            int tot_sample_id = 0;
            for(int mis_obj_idx=0; mis_obj_idx < mis_samples.size(); mis_obj_idx++) {
                for(int cur_sid=0; cur_sid < mis_samples[mis_obj_idx]; cur_sid++) {
                    int s_j = tot_sample_id / sample_len;
                    int s_i = tot_sample_id % sample_len;
                    // Sample a ray: stratified
                    double pixel_i = (s_i + random_double()) / sample_len - 0.5;
                    double pixel_j = (s_j + random_double()) / sample_len - 0.5;
                    // double pixel_i = random_double() - 0.5;
                    // double pixel_j = random_double() - 0.5;
                    Vec pixel_loc = pixel_orig + pixel_u * (i + pixel_i) + pixel_v * (j + pixel_j);
                    Ray sampled_ray(lookfrom, pixel_loc - lookfrom, random_double());
                    // Raytracing with BVH
                    Ray r_in = sampled_ray;
                    pixel_color = pixel_color + radiance(point_lights, objects, bvh_root, r_in, ray_tracing_depth, background_color, rl_init, mis_obj_idx, sample_len * sample_len, mis_samples);
                    tot_sample_id += 1;
                }
            }
            image[j * image_width + i] = get_color(pixel_color / double(sample_len * sample_len));
        }
    }
    // Store the image
    for (int i=0; i < image.size(); i++) {
        fprintf(out, "%d %d %d\n", int(image[i].x), int(image[i].y), int(image[i].z));
    }
    std::clog << "\rDone.                 \n";
}