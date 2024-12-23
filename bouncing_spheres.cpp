#include "smallpt.hpp"

// Source: https://raytracing.github.io/
// This is to test BVH optimization
void bouncing_spheres() {
    FILE *f;
    f = fopen("images/bouncing_spheres.ppm", "w");

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

    auto material1 = std::make_shared<TransmissiveMaterial>(1.5);
    world.push_back(std::make_shared<Sphere>(Vec(0, 1, 0), 1.0, material1));

    auto material2 = std::make_shared<DiffusiveMaterial>(Vec(0.4, 0.2, 0.1));
    world.push_back(std::make_shared<Sphere>(Vec(-4, 1, 0), 1.0, material2));

    auto material3 = std::make_shared<SpecularMaterial>(Vec(0.7, 0.6, 0.5), 0.0);
    world.push_back(std::make_shared<Sphere>(Vec(4, 1, 0), 1.0, material3));

    int image_width = 400;
    int image_height = 400 * 9 / 16;
    int sample_len = 10; // samples_per_pixel = sample_len * sample_len
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

int main() {
    bouncing_spheres();
}

// Reproduction: 
// $ g++ -O3 -fopenmp bouncing_spheres.cpp -o bouncing_spheres
// $ time ./bouncing_spheres