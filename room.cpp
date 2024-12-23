#include "smallpt.hpp"

// The custom scene emulating https://benedikt-bitterli.me/resources/ "Veach, Bidir Room by Benedikt Bitterli"
void room() {
    FILE *f;
    f = fopen("images/room.ppm", "w");
    std::vector<std::shared_ptr<Object>> world;
    std::vector<std::shared_ptr<Object>> lights;

    // Texture
    // Vec brown(0.8, 0.4, 0);
    Vec brown(0.0, 0.0, 0.0);
    Vec gray(0.5, 0.5, 0.5);
    Vec white(0.87, 0.87, 0.87);
    Vec light(1, 0.87, 0.6);
    auto empty_material = std::shared_ptr<Material>();

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
    lights.push_back(std::make_shared<Sphere>(toplight_orig, toplight_radius, empty_material));
    

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