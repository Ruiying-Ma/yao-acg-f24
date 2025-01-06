#include "smallpt.hpp"

// Source: https://benedikt-bitterli.me/resources/images/bedroom.png
void bedroom() {
    FILE *f;
    f = fopen("images/bedroom.ppm", "w");
    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;

    auto boxes = std::make_shared<DiffusiveMaterial>(Vec(0.48304399847984314, 0.38466399908065796, 0.30156099796295166));
    auto plastic_cable = std::make_shared<DiffusiveMaterial>(Vec(0.558543, 0.558543, 0.558543));
    auto lamp_emitter = std::make_shared<DiffusiveMaterial>(Vec(0.64, 0.64, 0.64));
    auto blankets = std::make_shared<DiffusiveMaterial>(Vec(0.48543548583984375, 0.4562634229660034, 0.42807474732398987));
    auto bed_sheets = std::make_shared<DiffusiveMaterial>("texture/bedroom/wallpaper-1.jpg");
    auto window = std::make_shared<DiffusiveMaterial>(Vec(0.028173, 0.028173, 0.028173));
    auto picture_backing = std::make_shared<DiffusiveMaterial>(Vec(0.11156699806451797, 0.037067998200654984, 0.01701599918305874));
    auto picture = std::make_shared<DiffusiveMaterial>("texture/bedroom/Teapot.png");
    auto rocks1 = std::make_shared<DiffusiveMaterial>(Vec(0.35082700848579407,0.24298599362373352,0.17882999777793884));
    auto rocks2 = std::make_shared<DiffusiveMaterial>(Vec(0.098964, 0.098964, 0.098964));
    auto rocks3 = std::make_shared<DiffusiveMaterial>(Vec(0.558544, 0.558544, 0.558544));
    auto deco_plant = std::make_shared<DiffusiveMaterial>(Vec(0.04177200049161911, 0.01130599994212389, 0.007575000170618296));
    auto painting = std::make_shared<DiffusiveMaterial>(Vec(0.015396, 0.015396, 0.015396));
    auto carpet = std::make_shared<DiffusiveMaterial>(Vec(0.034499, 0.034499, 0.034499));
    auto matress = std::make_shared<DiffusiveMaterial>(Vec(0.893289, 0.893289, 0.893289));
    auto wood_floor = std::make_shared<DiffusiveMaterial>("texture/bedroom/wood4.jpg");
    auto walls2 = std::make_shared<DiffusiveMaterial>(Vec(0.799999, 0.799999, 0.799999));
    auto wood_furniture = std::make_shared<DiffusiveMaterial>("texture/bedroom/panel-wood-3.jpg");
    auto walls = std::make_shared<DiffusiveMaterial>(Vec(0.799999, 0.799999, 0.799999));
    auto mirror = std::make_shared<SpecularMaterial>(Vec(1.0, 1.0, 1.0), 0.1);
    auto aluminum = std::make_shared<SpecularMaterial>(Vec(0.8, 0.8, 0.8), 0.2);
    auto book_cover = std::make_shared<DiffusiveMaterial>(Vec(0.0, 0.0, 0.0));
    auto book_pages = std::make_shared<DiffusiveMaterial>(Vec(0.567027, 0.567027, 0.567027));
    auto lamp_metal = std::make_shared<SpecularMaterial>(Vec(0.8, 0.8, 0.8), 0.1);
    auto lamp_glass = std::make_shared<TransmissiveMaterial>(1.5);
    auto picture_frame = std::make_shared<SpecularMaterial>(Vec(0.8, 0.8, 0.8), 0.1);
    auto glass = std::make_shared<TransmissiveMaterial>(1.5);
    auto vase = std::make_shared<TransmissiveMaterial>(1.5);
    auto curtains = std::make_shared<DiffusiveMaterial>(Vec(0.531049, 0.531049, 0.531049));
    auto curtain_rod = std::make_shared<SpecularMaterial>(Vec(0.5, 0.5, 0.5), 0.1);
    auto stainless_smooth = std::make_shared<SpecularMaterial>(Vec(1.0, 1.0, 1.0), 0.1);

    std::string relative_mesh_dir = "mesh/bedroom/";

    std::vector<std::vector<Tungsten::Vertex>> vtabs(68);
    std::vector<std::vector<Tungsten::TriangleI>> ftabs(68);

    std::vector<std::shared_ptr<Material>> mats = {
        // 0-9
        blankets, blankets, boxes, picture, curtains, stainless_smooth, wood_furniture, wood_furniture, walls2, wood_furniture,
        // 10-19
        window, mirror, wood_furniture, walls, curtain_rod, curtain_rod, blankets, painting, painting, carpet, 
        // 20-29
        carpet, matress, aluminum, lamp_metal, stainless_smooth, rocks1, lamp_glass, aluminum, stainless_smooth, stainless_smooth,
        // 30-39
        wood_furniture, wood_furniture, wood_furniture, deco_plant, bed_sheets, rocks3, aluminum, lamp_glass, mirror, window,
        // 40-49
        aluminum, bed_sheets, aluminum, aluminum, aluminum, stainless_smooth, aluminum, aluminum, lamp_glass, lamp_glass,
        // 50-59
        picture_frame, lamp_glass, glass, picture_backing, lamp_emitter, rocks2, lamp_metal, curtains, plastic_cable, lamp_emitter,
        // 60-68
        wood_floor, lamp_emitter, plastic_cable, book_cover, book_pages, vase, lamp_metal, vase, bed_sheets, 
    };
    assert(mats.size() == 69);

    for(int i = 0; i < 68; i++) {
        if (i == 4 || i == 57) {continue;}
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
        std::clog<<"Create mesh "<<mesh_id<<"\n";/////////////
        Vec trans = Vec();
        if (mesh_id == "20" || mesh_id == "19") {trans = Vec(-0.4418520927429199, 0.0, 0.0);}
        else if (mesh_id == "18") {trans = Vec( -0.10999999940395355, 0.0, 0.0);}
        else if (mesh_id == "17") {trans = Vec(-0.2358992099761963, 0.0, 0.0);}
        create_mesh(world, relative_mesh_path, vtabs[i], ftabs[i], mats[i], trans);
        if (mesh_id == "18") {
            trans = Vec(-1.1020985841751099, 0.09580427408218384, 0.0);
            create_mesh(world, relative_mesh_path, vtabs[i], ftabs[i], mats[i], trans);
        }
    }

    

    int image_width = 1280; //1280
    int image_height = 720; //720
    int sample_len = 28;
    Vec background = Vec(0.70, 0.80, 1.00);

    double vfov     = 65.0;
    Vec lookfrom = Vec(3.4555792808532715, 1.2124358415603638, 3.2989654541015625);
    Vec lookat   = Vec(0.0942695364356041, 1.1369876861572266, 0.39623117446899414);
    Vec vup      = Vec(0, 1, 0);
    double focus = 10;

    auto empty_material = std::shared_ptr<Material>();
    auto light = std::make_shared<AreaLight>(Vec(1.0, 1.0, 1.0));
    world.push_back(std::make_shared<Quad>(Vec(1.4437944889068604, 0.0,-1.2673544883728027), Vec(1.815588116645813 * 2, 0.0, 0.0), Vec(0.0, 1.0648202896118164 * 2.5, 0.0), light));
    world.push_back(std::make_shared<Quad>(Vec(1.4437944889068604, 0.0,-1.2673544883728027), Vec(0.0, 1.0648202896118164 * 2.5, 0.0), Vec(-1.815588116645813 * 2, 0.0, 0.0), light));
    lights.push_back(std::make_shared<Quad>(Vec(1.4437944889068604, 0.0,-1.2673544883728027), Vec(1.815588116645813 * 2, 0.0, 0.0), Vec(0.0, 1.0648202896118164 * 2.5, 0.0), empty_material));
    lights.push_back(std::make_shared<Quad>(Vec(1.4437944889068604, 0.0,-1.2673544883728027), Vec(0.0, 1.0648202896118164 * 2.5, 0.0), Vec(-1.815588116645813 * 2, 0.0, 0.0), empty_material));
    

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
    bedroom();
}

// Reproduction: 
// $ g++ -O3 -fopenmp bedroom.cpp -o bedroom
// $ time ./bedroom
