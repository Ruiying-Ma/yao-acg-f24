#ifndef CAMERA_H
#define CAMERA_H
#include "hittable.h"

class Camera {
  public:
    double aspect_ratio      = 1.0;  // Ratio of image width over height
    int    image_width       = 100;  // Rendered image width in pixel count
    int    samples_per_pixel = 10;   // Count of random samples for each pixel
    int    max_depth         = 10;   // Maximum number of ray bounces into scene
    color_t  background;             // Scene background color

    double vfov         = 90;                   // Vertical view angle (field of view)
    point_t lookfrom    = point_t(0,0,0);       // Point camera is looking from, equals camera.center_
    point_t lookat      = point_t(0,0,-1);      // Point camera is looking at
    direction_t   vup   = direction_t(0,1,0);   // Camera-relative "up" direction

    double defocus_angle = 0;  // Variation angle of rays through each pixel
    double focus_dist = 10;    // Distance from camera lookfrom point to plane of perfect focus

    void render(FILE *out, const Hittable& world) {
        initialize_();
        fprintf(out, "P3\n%d %d\n%d\n", image_width, image_height_, 255);
        std::vector<color_t> store_img(image_height_ * image_width);
        // color_t pixel_color(0,0,0);
        // #pragma omp parallel for schedule(dynamic, 1) private(pixel_color)
        for (int j = 0; j < image_height_; j++) {
            std::clog << "\rScanlines remaining: " << (image_height_ - j) << ' ' << std::flush;
            for (int i = 0; i < image_width; i++) {
                color_t pixel_color(0,0,0);
                // Anti-aliasing: uniform sampling
                for (int sample = 0; sample < samples_per_pixel; sample++) {
                    Ray r = get_ray_(i, j);
                    pixel_color += ray_color_(r, max_depth, world);
                }
                store_img[j * image_width + i] = get_color(pixel_samples_scale_ * pixel_color);
            }
            
        }

        // Store the image
        for (int i=0; i < store_img.size(); i++) {
            fprintf(out, "%d %d %d\n", int(store_img[i].x()), int(store_img[i].y()), int(store_img[i].z()));
        }

        std::clog << "\rDone.                 \n";
    }

  private:
    int    image_height_;           // Rendered image height
    double pixel_samples_scale_;    // Color scale factor for a sum of pixel samples; uniform sampling
    point_t center_;                // Camera center_
    point_t pixel00_loc_;           // Location of pixel 0, 0
    Vec   pixel_delta_u_;           // Offset to pixel to the right
    Vec   pixel_delta_v_;           // Offset to pixel below
    Vec   u_, v_, w_;               // Camera frame basis vectors; unit_length
    Vec   defocus_disk_u_;          // Defocus disk horizontal radius
    Vec   defocus_disk_v_;          // Defocus disk vertical radius

    void initialize_() {
        image_height_ = int(image_width / aspect_ratio);
        image_height_ = (image_height_ < 1) ? 1 : image_height_;

        pixel_samples_scale_ = 1.0 / samples_per_pixel;

        center_ = lookfrom;

        // Determine viewport dimensions.
        // viewport is at the perfect focus determined by focus_dist
        auto theta = degrees_to_radians(vfov);
        auto h = std::tan(theta/2);
        auto viewport_height = 2 * h * focus_dist;
        auto viewport_width = viewport_height * (double(image_width)/image_height_);

        // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
        w_ = unit_vector(lookfrom - lookat);
        u_ = unit_vector(cross(vup, w_));
        v_ = cross(w_, u_);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        Vec viewport_u = viewport_width * u_;    // Vector across viewport horizontal edge
        Vec viewport_v = viewport_height * -v_;  // Vector down viewport vertical edge

        // Calculate the horizontal and vertical delta vectors from pixel to pixel.
        pixel_delta_u_ = viewport_u / image_width;
        pixel_delta_v_ = viewport_v / image_height_;

        // Calculate the location of the upper left pixel.
        auto viewport_upper_left = center_ - (focus_dist * w_) - viewport_u/2 - viewport_v/2;
        pixel00_loc_ = viewport_upper_left + 0.5 * (pixel_delta_u_ + pixel_delta_v_);

        // Calculate the camera defocus disk basis vectors.
        auto defocus_radius = focus_dist * std::tan(degrees_to_radians(defocus_angle / 2));
        defocus_disk_u_ = u_ * defocus_radius;
        defocus_disk_v_ = v_ * defocus_radius;
    }

    // Construct a camera ray originating from the defocus disk and directed at a randomly sampled point around the pixel location i, j.
    Ray get_ray_(int i, int j) const {
        

        point_t offset = sample_point_from_square_();
        point_t pixel_sample = pixel00_loc_
                          + ((i + offset.x()) * pixel_delta_u_)
                          + ((j + offset.y()) * pixel_delta_v_);

        point_t ray_origin = (defocus_angle <= 0) ? center_ : sample_point_from_defocus_disk_();
        direction_t ray_direction = pixel_sample - ray_origin;
        double ray_time = random_double();

        return Ray(ray_origin, ray_direction, ray_time);
    }

    // Returns a random point in the [-.5,-.5]-[+.5,+.5] unit square.
    point_t sample_point_from_square_() const {
        return point_t(random_double() - 0.5, random_double() - 0.5, 0);
    }

    // Returns a random point in the disk of `radius` centered at the origin.
    point_t sample_point_from_disk_(double radius) const {
        return radius * random_in_unit_disk();
    }

    // Returns a random point inside the camera's defocus disk.
    point_t sample_point_from_defocus_disk_() const {
        point_t p = random_in_unit_disk();
        return center_ + (p[0] * defocus_disk_u_) + (p[1] * defocus_disk_v_);
    }

    // Return the color of the ray after scattered by the hittables in the world
    color_t ray_color_(const Ray& r, int depth, const Hittable& world) const {
        // If we've exceeded the ray bounce limit, no more light is gathered.
        if (depth <= 0)
            return color_t(0,0,0);

        HitRecord hit_rec;

        // If the ray hits nothing, return the background color.
        if (!world.hit(r, Interval(0.001, infinity), hit_rec))
            return background;

        Ray scattered;
        color_t attenuation;
        color_t color_from_emission = hit_rec.mat->emit(hit_rec.u, hit_rec.v, hit_rec.ip);

        // If the hittable doesn't scatter the ray
        if (!hit_rec.mat->scatter(r, hit_rec, attenuation, scattered))
            return color_from_emission;

        // Update the color by the attenuation from scattering
        color_t color_from_scatter = attenuation * ray_color_(scattered, depth-1, world);

        return color_from_emission + color_from_scatter;
    }
};


#endif
