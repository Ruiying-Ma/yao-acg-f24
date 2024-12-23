#include "smallpt.hpp"

void rocky_ground(bool use_normal_map, bool use_height_map) {
    FILE *f;
    std::string img_path = "images/rocky_ground";
    if (use_normal_map) {img_path += "_normal";}
    if (use_height_map) {img_path += "_height";}
    img_path += ".ppm";
    f = fopen(img_path.c_str(), "w");

    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;

    std::shared_ptr<DiffusiveMaterial> rocky_ground;
    if ((!use_normal_map) && (!use_height_map)) {rocky_ground = std::make_shared<DiffusiveMaterial>("texture/rocky_ground/color.jpg");}
    else if ((use_normal_map) && (!use_height_map)) {rocky_ground = std::make_shared<DiffusiveMaterial>("texture/rocky_ground/color.jpg", "texture/rocky_ground/normal.jpg");}
    else {{rocky_ground = std::make_shared<DiffusiveMaterial>("texture/rocky_ground/color.jpg", "texture/rocky_ground/normal.jpg", "texture/rocky_ground/height.jpg");}}

    auto light = std::make_shared<AreaLight>(Vec(1, 1, 1));
    auto empty_material = std::shared_ptr<Material>();

    world.push_back(std::make_shared<Quad>(Vec(10, 0, 10), Vec(0, 0, -20), Vec(-20, 20, 0), rocky_ground));
    world.push_back(std::make_shared<Quad>(Vec(40, 25, 40), Vec(-35, 0, 0), Vec(0, 0, -40), light));
    lights.push_back(std::make_shared<Quad>(Vec(40, 25, 40), Vec(-35, 0, 0), Vec(0, 0, -40), empty_material));

    int image_width = 400;
    int image_height = 400;
    int sample_len = 100;
    Vec background = Vec();

    double vfov     = 30;
    Vec lookfrom = Vec(60, 10, 0);
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
    rocky_ground(false, false);
}

// Reproduction: 
// $ g++ -O3 -fopenmp rocky_ground.cpp -o rocky_ground
// $ time ./rocky_ground