#include "include/camera.h"
#include "include/bvh.h"
#include <iostream>

void cornell_box() {
    
    HittableVector world;
    
    auto red   = std::make_shared<DiffuseMaterial>(color_t(.65, .05, .05));
    auto white = std::make_shared<DiffuseMaterial>(color_t(.73, .73, .73));
    auto green = std::make_shared<DiffuseMaterial>(color_t(.12, .45, .15));
    auto light = std::make_shared<AreaLight>(color_t(15, 15, 15));

    world.push_back(std::make_shared<Quad>(point_t(555,0,0), Vec(0,555,0), Vec(0,0,555), green));
    world.push_back(std::make_shared<Quad>(point_t(0,0,0), Vec(0,555,0), Vec(0,0,555), red));
    world.push_back(std::make_shared<Quad>(point_t(343, 554, 332), Vec(-130,0,0), Vec(0,0,-105), light));
    world.push_back(std::make_shared<Quad>(point_t(0,0,0), Vec(555,0,0), Vec(0,0,555), white));
    world.push_back(std::make_shared<Quad>(point_t(555,555,555), Vec(-555,0,0), Vec(0,0,-555), white));
    world.push_back(std::make_shared<Quad>(point_t(0,0,555), Vec(555,0,0), Vec(0,555,0), white));

    std::shared_ptr<Hittable> box1 = std::make_shared<Box>(point_t(0,0,0), point_t(165,330,165), white);
    box1 = std::make_shared<RotateYHittable>(box1, 15);
    box1 = std::make_shared<TranslateHittable>(box1, Vec(265,0,295));
    world.push_back(box1);

    std::shared_ptr<Hittable> box2 = std::make_shared<Box>(point_t(0,0,0), point_t(165,165,165), white);
    box2 = std::make_shared<RotateYHittable>(box2, -18);
    box2 = std::make_shared<TranslateHittable>(box2, Vec(130,0,65));
    world.push_back(box2);

    Camera cam;

    cam.aspect_ratio      = 1.0;
    cam.image_width       = 400;
    cam.samples_per_pixel = 16;
    cam.max_depth         = 5;
    cam.background        = color_t(0,0,0);

    cam.vfov     = 40;
    cam.lookfrom = point_t(278, 278, -800);
    cam.lookat   = point_t(278, 278, 0);
    cam.vup      = direction_t(0,1,0);

    cam.defocus_angle = 0;

    
    FILE *f = fopen("images/cornell_box.ppm", "w");

    world = HittableVector(std::make_shared<BVHNode>(world));

    cam.render(f, world);
}

void bouncing_spheres() {
    HittableVector world;

    world.push_back(std::make_shared<Sphere>(point_t(0,-1000,0), 1000, std::make_shared<DiffuseMaterial>(color_t(.2, .3, .1))));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point_t center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

            if ((center - point_t(4, 0.2, 0)).length() > 0.9) {
                std::shared_ptr<Material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color_t::random() * color_t::random();
                    sphere_material = std::make_shared<DiffuseMaterial>(albedo);
                    auto center2 = center + Vec(0, random_double(0,.5), 0);
                    world.push_back(std::make_shared<Sphere>(center, center2, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color_t::random(0.5, 1);
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
    world.push_back(std::make_shared<Sphere>(point_t(0, 1, 0), 1.0, material1));

    auto material2 = std::make_shared<DiffuseMaterial>(color_t(0.4, 0.2, 0.1));
    world.push_back(std::make_shared<Sphere>(point_t(-4, 1, 0), 1.0, material2));

    auto material3 = std::make_shared<SpecularMaterial>(color_t(0.7, 0.6, 0.5), 0.0);
    world.push_back(std::make_shared<Sphere>(point_t(4, 1, 0), 1.0, material3));

    world = HittableVector(std::make_shared<BVHNode>(world));

    Camera cam;

    cam.aspect_ratio      = 16.0 / 9.0;
    cam.image_width       = 400;
    cam.samples_per_pixel = 100;
    cam.max_depth         = 50;
    cam.background        = color_t(0.70, 0.80, 1.00);

    cam.vfov     = 20;
    cam.lookfrom = point_t(13,2,3);
    cam.lookat   = point_t(0,0,0);
    cam.vup      = direction_t(0,1,0);

    cam.defocus_angle = 0.6;
    cam.focus_dist    = 10.0;

    FILE *f = fopen("images/bouncing_spheres.ppm", "w");
    cam.render(f, world);
}


// g++ -O3 -fopenmp main.cpp -o main
int main() {
    // cornell_box();
    bouncing_spheres();
}