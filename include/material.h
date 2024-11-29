#ifndef MATERIAL_H
#define MATERIAL_H

#include "texture.h"

class Material;
class HitRecord;

class HitRecord {
  public:
    point_t ip; // intersection point
    direction_t normal; // **unit-length** normal of the intersected surface; same orientation as ray.
    std::shared_ptr<Material> mat; // material of the intersected surface
    double dist; // the distance from the ray origin to the intersection point
    double u; // u-coordinate in the texture space
    double v; // v-coordinate in the texture space
    bool front_face; // whether the intersected ray and the outward_normal is on the same orientation.

    // Set the hit record normal vector.
    // outward_normal: unit length.
    // Used by Hittable::hit
    void set_face_normal(const Ray& r, const direction_t& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

class Material {
    public:
    virtual ~Material() = default;

    // Return the emitted color. Default as black (0, 0, 0), i.e., no emitted light.
    virtual color_t emit(double u, double v, const point_t& p) const {
        return color_t(0,0,0);
    }

    // Return whether the material scatters the light. Default as false (e.g., for lights)
    // r_in: ray that comes in
    // hit_rec: the record of the current hit; calculated by Hittable::hit.
    // attenuation: **store** the attenuated factor for the color of r_in after scattered
    // scattered: **store** the scattered ray
    virtual bool scatter(
        const Ray& r_in, const HitRecord& hit_rec, color_t& attenuation, Ray& scattered
    ) const {
        return false;
    }
};

class DiffuseMaterial : public Material {
    public:
    // Define a colored diffuse material
    // albedo: attenuation
    DiffuseMaterial(const color_t& albedo) : texture_(std::make_shared<ColorTexture>(albedo)) {}

    DiffuseMaterial(std::shared_ptr<Texture> tex) : texture_(tex) {}

    // Scatter the incoming ray to a random direction
    bool scatter(const Ray& r_in, const HitRecord& hit_rec, color_t& attenuation, Ray& scattered)
    const override {
        // A random direction for scattered ray
        direction_t scatter_dir = hit_rec.normal + random_unit_vector();

        // Catch degenerate scatter direction
        if (scatter_dir.near_zero())
            scatter_dir = hit_rec.normal;

        // Store attenuation, scattered
        scattered = Ray(hit_rec.ip, scatter_dir, r_in.time());
        attenuation = texture_->value(hit_rec.u, hit_rec.v, hit_rec.ip);
        return true;
    }

  private:
    std::shared_ptr<Texture> texture_;
};


class SpecularMaterial : public Material {
    public:
    // Define a colored specular material
    // albedo: attenuation
    SpecularMaterial(const color_t& albedo, double fuzz) : texture_(std::make_shared<ColorTexture>(albedo)), fuzz_(fuzz < 1 ? fuzz : 1) {}

    bool scatter(const Ray& r_in, const HitRecord& hit_rec, color_t& attenuation, Ray& scattered)
    const override {
        direction_t reflected = unit_vector(reflect(r_in.direction(), hit_rec.normal));
        reflected = reflected + (fuzz_ * random_unit_vector());
        scattered = Ray(hit_rec.ip, reflected, r_in.time());
        attenuation = texture_->value(hit_rec.u, hit_rec.v, hit_rec.ip);
        return (dot(scattered.direction(), hit_rec.normal) > 0);
    }

  private:
    std::shared_ptr<Texture> texture_;
    double fuzz_;
};

// Currently doesn't support texture
class TransmissiveMaterial : public Material {
    public:
    TransmissiveMaterial(double refraction_index) : refraction_index_(refraction_index) {}

    bool scatter(const Ray& r_in, const HitRecord& hit_rec, color_t& attenuation, Ray& scattered)
    const override {
        attenuation = color_t(1.0, 1.0, 1.0); // no attenuation in refracted ray
        double ri = hit_rec.front_face ? (1.0/refraction_index_) : refraction_index_;

        direction_t unit_direction = unit_vector(r_in.direction());
        // theta: angle between r_in and -normal
        double cos_theta = std::fmin(dot(-unit_direction, hit_rec.normal), 1.0);
        double sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);

        bool cannot_refract = ri * sin_theta > 1.0;
        direction_t refracted_direction;

        if (cannot_refract || reflectance_(cos_theta, ri) > random_double())
            refracted_direction = reflect(unit_direction, hit_rec.normal); // total reflection
        else
            refracted_direction = refract(unit_direction, hit_rec.normal, ri);

        scattered = Ray(hit_rec.ip, refracted_direction, r_in.time());
        return true;
    }

  private:
    // Refractive index in vacuum or air, or the ratio of the material's refractive index over the refractive index of the enclosing media
    double refraction_index_;

    // Use Schlick's approximation for reflectance.
    static double reflectance_(double cosine, double refraction_index) {
        double r0 = (1 - refraction_index) / (1 + refraction_index);
        r0 = r0*r0;
        return r0 + (1-r0)*std::pow((1 - cosine),5);
    }
};


class AreaLight : public Material {
    public:
    AreaLight(std::shared_ptr<Texture> texture) : texture_(texture) {}
    // Create a colored area light
    AreaLight(const color_t& color) : texture_(std::make_shared<ColorTexture>(color)) {}

    color_t emit(double u, double v, const point_t& p) const override {
        return texture_->value(u, v, p);
    }

  private:
    std::shared_ptr<Texture> texture_;
};


// class DiffuseRandomSamplingMaterial : public Material {
//     public:
//     DiffuseRandomSamplingMaterial();

//     bool scatter(const Ray& r_in, const HitRecord& hit_rec, color_t& attenuation, Ray& scattered)
//     const override {
//         scattered = Ray(hit_rec.ip, )
//     }

//     private:
//     std::shared_ptr<Texture> texture_;
    
//     // Probaility distribution for scattering
//     // Return the probability of case (x, lambda, r_in , r_out) = (hit_rec, r_in, scattered)
//     double prob_scattering_(const Ray& r_in, const HitRecord& hit_rec, const Ray& scattered) const {
//         double cos_theta = dot(hit_rec.normal, unit_vector(scattered.direction()));
//         return cos_theta < 0 ? 0 : cos_theta;
//     }

//     double 
// }

#endif