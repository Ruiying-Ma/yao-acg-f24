# Midterm Report
All codes are written in `smallpt.cpp`.
## Demo Reproduction Guide
1. In `smallpt.cpp`, choose the demo in the main function. For example, if you want to render the cornell box scene:
    ```cpp
    int main() {
        cornell_box();
    }
    ```

    If you want to render the bouncing spheres:
    ```cpp
    int main() {
        bouncing_spheres();
    }
    ```

2. Create folder `images/` under the root.

3. Run cli
    ```
    $ g++ -O3 -fopenmp smallpt.cpp -o smallpt
    $ time ./smallpt
    ```

    If you render the cornell box scene, `cornell_box1.ppm` will be saved under `images/`. If you render the bouncing spheres, `bouncing_spheres1.ppm` will be saved under `images/`.

The BVH optimization can be seen using the bouncing spheres. The rendering time is about 18 sec with BVH. 

## References
1. I learned the code framework from [the textbook](https://raytracing.github.io/). My code in `smallpt.cpp` reimplemented it in a non-superficial way, except for:
    - the `ImageTexture` that reads the texture from an image. The three files
        - `external/stb_image_write.h`
        - `external/stb_image.h`
        - `include/rtw_stb_image.h`

        are directly copied from the textbook.

    Exemplar difference between my code and the textbook: 
    - The implementation of the textbook separates diffusive (lambertian) material from the specular and the transmissive ones by using `srec.skip_pdf`. This is unnecessary, and my code treats them equally. 
    - The textbook treats BVHNodes as objects (hittables). In contrast, I treat them solely as data structures for optimization. The optimized searching algorithm is independently implemented in function `intersect()`, rather than in class method. 

2. My code was built by expanding the 99-line [smallpt](https://www.kevinbeason.com/smallpt/). 