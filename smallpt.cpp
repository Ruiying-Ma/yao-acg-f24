#include "include/rtw_stb_image.h"

#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2009
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <cassert>
#include <iostream>
#include <memory>
#include <limits>
#include <vector>
#include <algorithm>


const double INF = std::numeric_limits<double>::infinity(); // for double
const int ray_tracing_depth = 1e5; // set this to be large, meaning that we are using roussian roulette for optimization
const double rl_init = 0.95; // the initial roussian roulette weight
const double cost_trav = 0.125; // the cost of traversal of SAH
const double cost_isect = 1.0; // the cost of ray intersecting an object in SAH
const int n_bucket = 12; // the number of buckets in SAH

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

struct Material;
struct HitRecord {
    Vec ip, n; // intersection_point, normal (reverse direction as r_in; always normalized)
    double u, v; // coordinate at the texture space
    bool front_face; // determine whether the object emits light
    double dist; // ip = ray.orig + n * dist
    std::shared_ptr<Material> mat;

    // Set `n` such that `dot(n, r_in) >= 0`
    void set_norm(const Ray& r_in, const Vec& out_norm_u) {
        front_face = r_in.d.dot(out_norm_u) < 0;
        n = (front_face) ? out_norm_u : out_norm_u * (-1);
    }

    void print() {
        std::clog<<"HitRecord:\n\tintersection: ";
        ip.print();
    }
};

struct Texture {
    virtual ~Texture() = default;
    virtual Vec at(double u, double v, const Vec& p) const = 0;
};

struct ColorTexture : Texture {
    Vec color;
    ColorTexture(const Vec& c_): color(c_) {}

    // Return the color at (u, v)
    Vec at(double u, double v, const Vec& p) const override {
        return color;
    }
};

// Copied from the textbook
struct ImageTexture : Texture {
    rtw_image img;
    ImageTexture(const char* img_path): img(img_path) {}
    Vec at(double u, double v, const Vec& p) const override {
        if(img.height() <=0) {return Vec(0, 1, 1);}
        u = clamp(u, 0.0, 1.0);
        v = 1 - clamp(v, 0.0, 1.0);
        auto pixel = img.pixel_data(int(u * img.width()), int(v * img.height()));
        return Vec(pixel[0] / 255.0, pixel[1] / 255.0, pixel[2] / 255.0);
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
};

struct DiffusiveMaterial : Material {
    std::shared_ptr<Texture> tex;
    DiffusiveMaterial(const Vec& c_): tex(std::make_shared<ColorTexture>(c_)) {}
    DiffusiveMaterial(const char* img_path): tex(std::make_shared<ImageTexture>(img_path)) {}
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
};

struct SpecularMaterial : Material {
    std::shared_ptr<Texture> tex;
    double fuzz;
    SpecularMaterial(const Vec& c_): tex(std::make_shared<ColorTexture>(c_)), fuzz(0.5) {}
    SpecularMaterial(const char* img_path): tex(std::make_shared<ImageTexture>(img_path)) {}
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
};

struct HDRLight : Material {
    std::shared_ptr<Texture> tex;
    HDRLight(const Vec& c_): tex(std::make_shared<ColorTexture>(c_)) {}
    HDRLight(const char* img_path): tex(std::make_shared<ImageTexture>(img_path)) {}
    double pScatter(const Ray& r_in, const Ray& r_out, const HitRecord& hit_rec) const override {assert(false); return 0.0;}
    Vec emit(const HitRecord& hit_rec) const override {
        if (hit_rec.front_face){
            std::cerr<<"HDRLight cannot have front_face=True\n";
            exit(1);
            return Vec();
        }
        return tex->at(hit_rec.u, hit_rec.v, hit_rec.ip);
    }
    Vec attenuation(const HitRecord& hit_rec) const override {assert(false);return Vec(1.0, 1.0, 1.0);}  
    Vec sample(const Ray& r_in, const HitRecord& hit_rec) const override {return Vec();} 
    double prob(const HitRecord& hit_rec, const Vec& direction) const override {return 0.0;}
    double sample_thresh() const override {return 0.0;}
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

    // Static sphere
    Sphere(const Vec& static_center, double radius, std::shared_ptr<Material> material) : o(Ray(static_center, Vec())), rad(std::fmax(radius, 0)), mat(material) {
        Vec radius_v = Vec(rad, rad, rad);
        bbox = BBox(static_center - radius_v, static_center + radius_v);
        /////////////////////////
        // std::clog<<"Create a sphere at ";
        // static_center.print();
        // std::clog<<"\n\tbbox: ";
        // bbox.print();
        // std::clog<<"\n";
        ////////////////////////////
    }

    // Moving sphere
    Sphere(const Vec& from, const Vec& to, double radius, std::shared_ptr<Material> material) : o(Ray(from, to - from)), rad(std::fmax(radius, 0)), mat(material) {
        Vec radius_v = Vec(rad, rad, rad);
        bbox = BBox(BBox(from - radius_v, from + radius_v), BBox(to - radius_v, to + radius_v));
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
        hit_rec.u = phi / (2 * M_PI);
        hit_rec.v = theta / M_PI;
        hit_rec.dist = root;
        hit_rec.ip = ip;
        hit_rec.set_norm(r, out_norm_u);
        hit_rec.mat = mat;
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
        Vec src2sphere = o.at(0) - src;
        double cos_theta = std::sqrt(1 - rad * rad / (src2sphere.dot(src2sphere)));
        return 1 / (2 * M_PI * (1 - cos_theta));
    }

    // Importance sampling; sample a point on the shape
    // Return: a sampled `direction`
    // Only for static sphere
    Vec sample(const Vec& src) const override {
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
        /////////////////////////
        // std::clog<<"Create a quad at ";
        // o.print();
        // u.print();
        // v.print();
        // std::clog<<"\n\tbbox: ";
        // bbox.print();
        // std::clog<<"\n";
        ////////////////////////////
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
        hit_rec.dist = dist;
        hit_rec.ip = ip;
        hit_rec.set_norm(r, n_u);
        hit_rec.mat = mat;
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
            std::vector<BucketInfo> buckets(n_bucket);
            double left_margin = objects[start_idx]->get_bbox().centroid(max_axis);
            double right_margin = objects[end_idx - 1]->get_bbox().centroid(max_axis);
            double bucket_size = (right_margin - left_margin) / n_bucket;
            assert(bucket_size > 0);
            //Initialize the buckets: O(n)
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
        left = std::make_shared<BVHNode>(objects, start_idx, mid);
        right = std::make_shared<BVHNode>(objects, mid, end_idx);
    }
    // Sort according to centroid, due to SAH
    static bool box_x_compare (const std::shared_ptr<Object> a, const std::shared_ptr<Object> b) {return a->get_bbox().centroid(0) < b->get_bbox().centroid(0);}
    static bool box_y_compare (const std::shared_ptr<Object> a, const std::shared_ptr<Object> b) {return a->get_bbox().centroid(1) < b->get_bbox().centroid(1);}
    static bool box_z_compare (const std::shared_ptr<Object> a, const std::shared_ptr<Object> b) {return a->get_bbox().centroid(2) < b->get_bbox().centroid(2);}
};

// Importance sampling
// Return the probability of the sampled ray
// Store the sampled ray to `r_out`
double sample(const Ray& r_in, const HitRecord& hit_rec, const std::vector<std::shared_ptr<Object>>& point_lights, Ray& r_out) {
    // Sample a scattering direction
    r_out.o = hit_rec.ip;
    r_out.t = r_in.t;
    double thresh = hit_rec.mat->sample_thresh();
    if (point_lights.size() == 0) {
        thresh = 0;
    }
    if (random_double() < thresh) {
        // sample using lights based on their shapes: randomly choose a light
        int pl_idx = random_int(0, point_lights.size() - 1);
        r_out.d = point_lights[pl_idx]->sample(hit_rec.ip); // src = intersection
    }
    else {
        // sample using the material of the hit object
        r_out.d = hit_rec.mat->sample(r_in, hit_rec);
    }
    // Calculate the sampling probability
    double pl_prob = 0.0;
    if (thresh > 0) {
        for (const auto& pl : point_lights) {
            pl_prob += pl->prob(hit_rec.ip, r_out.d); // origin = intersection, direction = out_dir
        }
        pl_prob /= point_lights.size();
    }
    double mat_prob = 0.0;
    if (1 - thresh > 0) {
        mat_prob = hit_rec.mat->prob(hit_rec, r_out.d);
    }
    return thresh * pl_prob + (1 - thresh) * mat_prob;
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
    double rl_weight
) {
    if (depth <= 0) {return Vec();}
    HitRecord hit_rec;
    int obj_idx = intersect(objects, bvh_root, r_in, hit_rec, 1e-3, INF);
    if (obj_idx == -1) {return background;}
    // std::clog<<"intersect with "<<obj_idx<<std::endl;///////////////////////////////////////
    // Roussian roulette
    if (random_double() > rl_weight) {return Vec();}
    // Sample the scattered ray
    Ray r_out;
    double sampled_prob = sample(r_in, hit_rec, point_lights, r_out);
    /////////////////////////////////////////
    // if (sampled_prob == 0) {
    //     std::clog<<"emit "; 
    //     hit_rec.mat->emit(hit_rec).print(); 
    //     std::clog<<std::endl;
    // }
    //////////////////////////////////////////
    if (sampled_prob == 0) {return hit_rec.mat->emit(hit_rec);}
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
    Vec post_color = radiance(point_lights, objects, bvh_root, r_out, depth - 1, background, new_rl_weight);
    // Calculate the final color
    return (emit + attenuation.mult(post_color) * scattered_prob / sampled_prob) / rl_weight;
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
    for(int j = 0; j < image_height; j++) {
        std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
        for(int i = 0; i < image_width; i++) {
            Vec pixel_color;
            // Anti-aliasing by sampling rays for each pixel
            for(int s_j = 0; s_j < sample_len; s_j++) {
                for(int s_i = 0; s_i < sample_len; s_i++) {
                    // Sample a ray: stratified
                    double pixel_i = (s_i + random_double()) / sample_len - 0.5;
                    double pixel_j = (s_j + random_double()) / sample_len - 0.5;
                    // double pixel_i = random_double() - 0.5;
                    // double pixel_j = random_double() - 0.5;
                    Vec pixel_loc = pixel_orig + pixel_u * (i + pixel_i) + pixel_v * (j + pixel_j);
                    Ray sampled_ray(lookfrom, pixel_loc - lookfrom, random_double());
                    // Raytracing with BVH
                    Ray r_in = sampled_ray;
                    pixel_color = pixel_color + radiance(point_lights, objects, bvh_root, r_in, ray_tracing_depth, background_color, rl_init);
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


void demo(bool has_point_light) {
    FILE *f;
    if (has_point_light) {f = fopen("images/demo1_is.ppm", "w");} 
    else {f = fopen("images/demo1_nois.ppm", "w");}

    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;

    // auto red   = std::make_shared<SpecularMaterial>(Vec(.65, .05, .05));
    auto red   = std::make_shared<TransmissiveMaterial>(0.5);
    auto white = std::make_shared<DiffusiveMaterial>(Vec(.73, .73, .73));
    auto green = std::make_shared<DiffusiveMaterial>(Vec(.12, .45, .15));
    auto light = std::make_shared<AreaLight>(Vec(15, 15, 15));
    auto empty_material = std::shared_ptr<Material>();
    
    world.push_back(std::make_shared<Sphere>(Vec(0, 10, 0), 10, red));
    world.push_back(std::make_shared<Sphere>(Vec(0, -1e4, 0), 1e4, green));
    world.push_back(std::make_shared<Sphere>(Vec(4, 20, 14), 2, light));
    if (has_point_light) {
        lights.push_back(std::make_shared<Sphere>(Vec(4, 20, 14), 2, empty_material));
    }

    int image_width = 400;
    int image_height = 400;
    int sample_len = 40;
    Vec background = Vec();

    double vfov     = 40;
    Vec lookfrom = Vec(60, 10, 20);
    Vec lookat   = Vec(0, 10, 0);
    Vec vup      = Vec(0, 1, 0);
    double focus = 10;


    render(
        f,
        lights, 
        world,
        image_height,
        image_width,
        sample_len,
        lookfrom, // camera position
        lookat, // camera target
        vup, // up direction for camera
        vfov, // camera vfov
        focus, // focus of the camera
        background
    );
}

void demo_hdr() {
    FILE *f;
    f = fopen("images/demo_hdr1.ppm", "w");

    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;

    auto red   = std::make_shared<SpecularMaterial>(Vec(0.85, 0.85, 0.85), 0.0);
    auto white = std::make_shared<DiffusiveMaterial>(Vec(.73, .73, .73));
    auto green = std::make_shared<DiffusiveMaterial>(Vec(.12, .45, .15));
    // auto light = std::make_shared<AreaLight>("texture/forest.jpg");
    auto light = std::make_shared<HDRLight>("texture/skysphere.jpg");
    auto empty_material = std::shared_ptr<Material>();
    
    world.push_back(std::make_shared<Sphere>(Vec(0, 10, 0), 10, red));
    // world.push_back(std::make_shared<Sphere>(Vec(0, -1e4, 0), 1e4, green));
    world.push_back(std::make_shared<Sphere>(Vec(0, 10, 0), 1e4, light));

    int image_width = 400;
    int image_height = 400;
    int sample_len = 10;
    Vec background = Vec();

    double vfov     = 40;
    Vec lookfrom = Vec(0, 10, 60);
    Vec lookat   = Vec(0, 10, 0);
    Vec vup      = Vec(0, 1, 0);
    double focus = 10;


    render(
        f,
        lights, 
        world,
        image_height,
        image_width,
        sample_len,
        lookfrom, // camera position
        lookat, // camera target
        vup, // up direction for camera
        vfov, // camera vfov
        focus, // focus of the camera
        background
    );
}

// This is to test BVH optimization
void bouncing_spheres() {
    FILE *f;
    f = fopen("images/bouncing_spheres1.ppm", "w");

    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;

    world.push_back(std::make_shared<Sphere>(Vec(0,-1000,0), 1000, std::make_shared<DiffusiveMaterial>(Vec(.2, .3, .1))));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            Vec center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

            if ((center - Vec(4, 0.2, 0)).len() > 0.9) {
                std::shared_ptr<Material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = Vec(random_double(), random_double(), random_double()).mult(Vec(random_double(), random_double(), random_double()));
                    sphere_material = std::make_shared<DiffusiveMaterial>(albedo);
                    auto center2 = center + Vec(0, random_double(0,.5), 0);
                    world.push_back(std::make_shared<Sphere>(center, center2, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = Vec(random_double(), random_double(), random_double());
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = std::make_shared<SpecularMaterial>(albedo, fuzz);
                    world.push_back(std::make_shared<Sphere>(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = std::make_shared<TransmissiveMaterial>(1.5);
                    world.push_back(std::make_shared<Sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = std::make_shared<TransmissiveMaterial>(0.5);
    world.push_back(std::make_shared<Sphere>(Vec(0, 1, 0), 1.0, material1));

    auto material2 = std::make_shared<DiffusiveMaterial>(Vec(0.4, 0.2, 0.1));
    world.push_back(std::make_shared<Sphere>(Vec(-4, 1, 0), 1.0, material2));

    auto material3 = std::make_shared<SpecularMaterial>(Vec(0.7, 0.6, 0.5), 0.0);
    world.push_back(std::make_shared<Sphere>(Vec(4, 1, 0), 1.0, material3));

    int image_width = 400;
    int image_height = 400 * 9 / 16;
    int sample_len = 10;
    Vec background = Vec(0.70, 0.80, 1.00);

    double vfov     = 20;
    Vec lookfrom = Vec(13,2,3);
    Vec lookat   = Vec(0, 0, 0);
    Vec vup      = Vec(0, 1, 0);
    double focus = 10;

    render(
        f,
        lights, 
        world,
        image_height,
        image_width,
        sample_len,
        lookfrom, // camera position
        lookat, // camera target
        vup, // up direction for camera
        vfov, // camera vfov
        focus, // focus of the camera
        background
    );
}

void cornell_box(bool has_point_light) {
    FILE *f;
    if (has_point_light) {f = fopen("images/cornell_box1_is.ppm", "w");} 
    else {f = fopen("images/cornell_box1_nois.ppm", "w");}

    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;
    
    auto red   = std::make_shared<SpecularMaterial>(Vec(.65, .05, .05));
    auto white = std::make_shared<DiffusiveMaterial>(Vec(.73, .73, .73));
    auto green = std::make_shared<DiffusiveMaterial>(Vec(.12, .45, .15));
    auto light = std::make_shared<AreaLight>(Vec(15, 15, 15));
    auto empty_material = std::shared_ptr<Material>();

    world.push_back(std::make_shared<Quad>(Vec(555,0,0), Vec(0,555,0), Vec(0,0,555), green));
    world.push_back(std::make_shared<Quad>(Vec(0,0,0), Vec(0,555,0), Vec(0,0,555), red));
    world.push_back(std::make_shared<Quad>(Vec(343, 554, 332), Vec(-130,0,0), Vec(0,0,-105), light));
    if (has_point_light) {
        lights.push_back(std::make_shared<Quad>(Vec(343, 554, 332), Vec(-130,0,0), Vec(0,0,-105), empty_material));
    }
    world.push_back(std::make_shared<Quad>(Vec(0,0,0), Vec(555,0,0), Vec(0,0,555), white));
    world.push_back(std::make_shared<Quad>(Vec(555,555,555), Vec(-555,0,0), Vec(0,0,-555), white));
    world.push_back(std::make_shared<Quad>(Vec(0,0,555), Vec(555,0,0), Vec(0,555,0), white));

    create_box(world, Vec(265, 0, 295), 165, 330, 165, white);
    
    create_box(world, Vec(130, 0, 65), 165, 165, 165, white);
    

    int image_width = 400;
    int image_height = 400;
    int sample_len = 10;
    Vec background = Vec();

    double vfov     = 40;
    Vec lookfrom = Vec(278, 278, -800);
    Vec lookat   = Vec(278, 278, 0);
    Vec vup      = Vec(0, 1, 0);
    double focus = 10;


    render(
        f,
        lights, 
        world,
        image_height,
        image_width,
        sample_len,
        lookfrom, // camera position
        lookat, // camera target
        vup, // up direction for camera
        vfov, // camera vfov
        focus, // focus of the camera
        background
    );
}


// The custom scene emulating https://benedikt-bitterli.me/resources/ "Veach, Bidir Room by Benedikt Bitterli"
void custom() {
    FILE *f;
    f = fopen("images/custom1.ppm", "w");
    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;

    // Texture
    // Vec brown(0.8, 0.4, 0);
    Vec brown(0.0, 0.0, 0.0);
    Vec gray(0.5, 0.5, 0.5);
    Vec white(0.87, 0.87, 0.87);
    Vec light(1, 0.87, 0.6);

    // The room: box, diffusive, white
    Vec room_orig(0, 0, 0);
    double room_y = 20;
    double room_z = room_y * 2.0;
    double room_x = room_z * 1.5;
    create_box(world, room_orig, room_x, room_y, room_z, std::make_shared<DiffusiveMaterial>(white));
    // The carpet: box, diffusive, carpet.jpg
    double carpet_orig_off_delta = room_z * 0.05;
    Vec carpet_orig = room_orig + Vec(carpet_orig_off_delta, 0.0, carpet_orig_off_delta);
    double carpet_z = room_z - 2.0 * carpet_orig_off_delta; // carpet thickness
    double carpet_x = carpet_z / 2.0;
    double carpet_y = 0.1;
    create_box(world, carpet_orig, carpet_x, carpet_y, carpet_z, std::make_shared<DiffusiveMaterial>("texture/carpet.jpg"));
    // The desk: boxes, diffusive, wood.jpg
    double desk_thickness = 0.4;
    double desk_ratio_x = 0.8;
    double desk_ratio_z = 0.6;
    double desk_face_off_y = room_y * 0.3; // desk height
    double desk_face_off_x = carpet_x * (1 - desk_ratio_x) * 0.5;
    double desk_face_off_z = carpet_z * (1 - desk_ratio_z) * 0.5;
    Vec desk_face_orig = carpet_orig + Vec(desk_face_off_x, desk_face_off_y, desk_face_off_z);
    create_box(world, desk_face_orig, carpet_x * desk_ratio_x, desk_thickness, carpet_z * desk_ratio_z, std::make_shared<DiffusiveMaterial>("texture/wood.jpg"));
    double desk_foot_off_x = desk_thickness * 1.5;
    double desk_foot_off_z = desk_face_off_x;
    create_box(world, desk_face_orig + Vec(desk_foot_off_x, -desk_face_off_y, desk_foot_off_z), desk_thickness, desk_face_off_y, desk_thickness, std::make_shared<DiffusiveMaterial>(gray));
    create_box(world, desk_face_orig + Vec(carpet_x * desk_ratio_x - desk_foot_off_x - desk_thickness, -desk_face_off_y, desk_foot_off_z), desk_thickness, desk_face_off_y, desk_thickness, std::make_shared<DiffusiveMaterial>(gray));
    create_box(world, desk_face_orig + Vec(desk_foot_off_x, -desk_face_off_y, carpet_z * desk_ratio_z - desk_foot_off_z - desk_thickness), desk_thickness, desk_face_off_y, desk_thickness, std::make_shared<DiffusiveMaterial>(gray));
    create_box(world, desk_face_orig + Vec(carpet_x * desk_ratio_x - desk_foot_off_x - desk_thickness, -desk_face_off_y, carpet_z * desk_ratio_z - desk_foot_off_z - desk_thickness), desk_thickness, desk_face_off_y, desk_thickness, std::make_shared<DiffusiveMaterial>(gray));
    // The football: sphere, diffusive, football.jpg
    double football_radius = desk_face_off_y * 0.25;
    Vec football_orig = room_orig + Vec(carpet_orig_off_delta + carpet_x * 1.2, football_radius, room_z * 0.75);
    world.push_back(std::make_shared<Sphere>(football_orig, football_radius, std::make_shared<DiffusiveMaterial>("texture/football.jpg")));
    // The frame: boxes, diffusive, wood.jpg
    double frame_remaining_y = room_y - carpet_y - desk_thickness - desk_face_off_y;
    double frame_orig_off_delta_y = frame_remaining_y * 0.2;
    double frame_y = frame_remaining_y - 2.0 * frame_orig_off_delta_y;
    double frame_z = frame_y * 16.0 / 9.0;
    double fram_orig_off_delta_z = (room_z - frame_z) / 2.0;
    Vec frame_orig = room_orig + Vec(0, carpet_y + desk_thickness + desk_face_off_y + frame_orig_off_delta_y, fram_orig_off_delta_z);
    double frame_thickness = desk_thickness * 0.5;
    create_box(world, frame_orig, frame_thickness, frame_y, frame_thickness, std::make_shared<DiffusiveMaterial>(brown));
    create_box(world, frame_orig + Vec(0, 0, frame_z - frame_thickness), frame_thickness, frame_y, frame_thickness, std::make_shared<DiffusiveMaterial>(brown));
    create_box(world, frame_orig + Vec(0, 0, frame_thickness), frame_thickness, frame_thickness, frame_z - 2 * frame_thickness, std::make_shared<DiffusiveMaterial>(brown));
    create_box(world, frame_orig + Vec(0, frame_y - frame_thickness, frame_thickness), frame_thickness, frame_thickness, frame_z - 2 * frame_thickness, std::make_shared<DiffusiveMaterial>(brown));
    // The photo: quad, diffusive, photo.jpg
    double photo_y = frame_y - 2 * frame_thickness;
    double photo_z = frame_z - 2 * frame_thickness;
    Vec photo_orig = frame_orig + Vec(frame_thickness * 0.2, frame_thickness, frame_thickness);
    world.push_back(std::make_shared<Quad>(photo_orig, Vec(0, photo_y, 0), Vec(0, 0, photo_z), std::make_shared<DiffusiveMaterial>("texture/photo.jpg")));
    // The top lights: sphere, area_light, light
    double toplight_radius = carpet_z * 0.5;
    Vec toplight_orig = room_orig + Vec(carpet_orig_off_delta + carpet_x * 1.1, room_y + toplight_radius * 0.9, room_z * 0.5);
    world.push_back(std::make_shared<Sphere>(toplight_orig, toplight_radius, std::make_shared<AreaLight>(light)));
    

    int image_width = 400;
    int image_height = 400 / 2;
    int sample_len = 100;
    Vec background = Vec(0.70, 0.80, 1.00);

    double vfov     = 40;
    Vec lookfrom = Vec(room_x * 0.95, room_y * 0.5, room_z * 0.5);
    Vec lookat   = Vec(room_x * 0.05, room_y * 0.5, room_z * 0.5);
    Vec vup      = Vec(0, 1, 0);
    double focus = 10;

    render(
        f,
        lights, 
        world,
        image_height,
        image_width,
        sample_len,
        lookfrom, // camera position
        lookat, // camera target
        vup, // up direction for camera
        vfov, // camera vfov
        focus, // focus of the camera
        background
    );
}

int main() {
    // bouncing_spheres();
    // demo(true);
    demo_hdr();
    // custom();
    // cornell_box(false);
}