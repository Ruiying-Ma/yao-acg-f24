#ifndef TEXTURE_H
#define TEXTURE_H

#include "ray.h"

class Texture {
    public:
    virtual ~Texture() = default;

    // Return the texture value at texture mapping (u, v)
    virtual texture_t value(double u, double v, const point_t& p) const = 0;
};


class ColorTexture : public Texture {
    public:
    ColorTexture(const color_t& rgb) : rgb_(rgb) {}

    ColorTexture(double r, double g, double b) : ColorTexture(color_t(r, g, b)) {}

    color_t value(double u, double v, const point_t& p) const override {
        return rgb_;
    }

    private:
    color_t rgb_;
};

#endif