#include "smallpt.hpp"

// Source: https://benedikt-bitterli.me/resources/images/car2.png
void car() {
    FILE *f;
    f = fopen("images/car.ppm", "w");
    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;

    auto car_paint = std::make_shared<SpecularMaterial>(Vec(.12, .45, .15), 0.1);
    auto window_glass = std::make_shared<TransmissiveMaterial>(1.5);
    auto black_rubber = std::make_shared<DiffusiveMaterial>(Vec(0, 0, 0));
    auto black = std::make_shared<DiffusiveMaterial>(Vec(0, 0, 0));
    auto white_rubber = std::make_shared<DiffusiveMaterial>(Vec(1, 1, 1));
    auto steel = std::make_shared<SpecularMaterial>(Vec(.8, .8, .8));
    auto ground = std::make_shared<DiffusiveMaterial>(Vec(.375, .375, .375));
    auto leather = std::make_shared<DiffusiveMaterial>(Vec(0.61, 0.36, 0.11));
    auto inner_body = std::make_shared<DiffusiveMaterial>(Vec(.25, .25, .25));
    auto dash = std::make_shared<DiffusiveMaterial>(Vec(.75, .75, .75));
    auto cabin = std::make_shared<DiffusiveMaterial>(Vec(0.15, 0.15, 0.15));
    auto chrome = std::make_shared<SpecularMaterial>(Vec(.5, .5, .5));
    auto hdr_light = std::make_shared<AreaLight>("texture/rosendal_park_sunset_puresky_8k.hdr");
    auto empty_material = std::shared_ptr<Material>();

    std::string relative_mesh_dir = "mesh/car/";

    std::vector<std::vector<Tungsten::Vertex>> vtabs(63);
    std::vector<std::vector<Tungsten::TriangleI>> ftabs(63);

    std::vector<std::shared_ptr<Material>> mats = {
        // 0-9
        steel, steel, steel, window_glass, window_glass, steel, steel, steel, steel, steel,
        // 10-19
        steel, steel, car_paint, white_rubber, black_rubber, steel, window_glass, window_glass, steel, inner_body, 
        // 20-29
        black, chrome, chrome, leather, chrome, white_rubber, steel, leather, chrome, black,
        // 30-39
        steel, chrome, black_rubber, steel, ground, chrome, steel, steel, steel, inner_body,
        // 40-49
        chrome, chrome, chrome, steel, car_paint, steel, car_paint, car_paint, inner_body, steel,
        // 50-59
        chrome, chrome, chrome, steel, steel, steel, steel, chrome, steel, dash,
        // 60-62
        car_paint, cabin, car_paint
    };

    for(int i = 0; i < 63; i++) {
        std::string mesh_id = std::to_string(i);
        int mesh_id_str_len = mesh_id.length();
        switch (mesh_id_str_len)
        {
        case 1:
            mesh_id = "00" + mesh_id;
            break;
        case 2:
            mesh_id = "0" + mesh_id;
            break;
        case 3:
            break;
        default:
            std::cerr<<"Invalid mesh id "<<mesh_id<<"\n";
            exit(1);
        }
        std::string relative_mesh_path = relative_mesh_dir + "Mesh" + mesh_id + ".wo3";
        create_mesh(world, relative_mesh_path, vtabs[i], ftabs[i], mats[i]);
    }

    

    int image_width = 1280; //1280
    int image_height = 720; //720
    int sample_len = 50;
    Vec background = Vec(0.70, 0.80, 1.00);

    double vfov     = 35.0;
    Vec lookfrom = Vec(-8.83707046508789, 5.837699890136719, 14.620699882507324);
    Vec lookat   = Vec(-1.8855600357055664, 1.7409499883651733, 2.2357900142669678);
    Vec vup      = Vec(0, 1, 0);
    double focus = 10;

    double hdr_rad = 1e4;
    world.push_back(std::make_shared<Sphere>(lookfrom, hdr_rad, hdr_light, true));
    lights.push_back(std::make_shared<Sphere>(lookfrom, hdr_rad, empty_material, true));

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
    car();
}

// Reproduction: 
// $ g++ -O3 -fopenmp car.cpp -o car
// $ time ./car
