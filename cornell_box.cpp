#include "smallpt.hpp"

void cornell_box(bool multiple_importance_sampling) {
    FILE *f;
    if (multiple_importance_sampling) {f = fopen("images/cornell_box_mis.ppm", "w");} 
    else {f = fopen("images/cornell_box_nomis.ppm", "w");}

    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;
    
    auto red   = std::make_shared<DiffusiveMaterial>(Vec(.65, .05, .05));
    auto white = std::make_shared<DiffusiveMaterial>(Vec(.73, .73, .73));
    auto green = std::make_shared<DiffusiveMaterial>(Vec(.12, .45, .15));
    auto light = std::make_shared<AreaLight>(Vec(15, 15, 15));
    auto empty_material = std::shared_ptr<Material>();

    world.push_back(std::make_shared<Quad>(Vec(555,0,0), Vec(0,555,0), Vec(0,0,555), green));
    world.push_back(std::make_shared<Quad>(Vec(0,0,0), Vec(0,555,0), Vec(0,0,555), red));
    world.push_back(std::make_shared<Quad>(Vec(343, 554, 332), Vec(-130,0,0), Vec(0,0,-105), light));
    if (multiple_importance_sampling) {
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

int main() {
    cornell_box(true);
}

// Reproduction: 
// $ g++ -O3 -fopenmp cornell_box.cpp -o cornell_box
// $ time ./cornell_box