#include "smallpt.hpp"

void transmissive_sphere(bool multiple_importance_sampling) {
    FILE *f;
    if (multiple_importance_sampling) {f = fopen("images/transmissive_sphere_mis.ppm", "w");} 
    else {f = fopen("images/transmissive_sphere_nomis.ppm", "w");}

    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;

    auto glass   = std::make_shared<TransmissiveMaterial>(0.5);
    auto green = std::make_shared<DiffusiveMaterial>(Vec(.12, .45, .15));
    auto light = std::make_shared<AreaLight>(Vec(15, 15, 15));
    auto empty_material = std::shared_ptr<Material>();
    
    world.push_back(std::make_shared<Sphere>(Vec(0, 10, 0), 10, glass));
    world.push_back(std::make_shared<Quad>(Vec(500, 0, 500), Vec(0, 0, -1000), Vec(-1000, 0, 0), green));
    world.push_back(std::make_shared<Sphere>(Vec(4, 20, 14), 2, light));
    if (multiple_importance_sampling) {
        lights.push_back(std::make_shared<Sphere>(Vec(4, 20, 14), 2, empty_material));
    }

    int image_width = 400;
    int image_height = 400;
    int sample_len = 10;
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

int main() {
    transmissive_sphere(false);
}

// Reproduction: 
// $ g++ -O3 -fopenmp transmissive_sphere.cpp -o transmissive_sphere
// $ time ./transmissive_sphere