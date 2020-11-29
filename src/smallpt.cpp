#include <bits/stdc++.h>

#include "FloatUtil.hpp"

#include "Timer.hpp"
#include "Random.hpp"
#include "Vector.hpp"
#include "Ray.hpp"
#include "Material.hpp"
#include "Sphere.hpp"
#include "Scene.hpp"
#include "PathTracing.hpp"

int main(int argc, char *argv[])
{
    srand(time(NULL));

    Timer timer_progress;
    Timer timer_rendertime;

    int img_width = 256;
    int img_height = 192;
    int sample_per_subpixel = 8;                                          // # samples
    Ray camera(Vector(50, 52, 295.6), Vector(0, -0.042612, -1).Normal()); // cam pos, dir
    Vector film_x_axis = Vector(img_width * .5135 / img_height);
    Vector film_y_axis = (film_x_axis % camera.direction).Normal() * .5135;
    Vector *img_buffer = new Vector[img_width * img_height];
#pragma omp parallel for schedule(dynamic, 1) private(ray) // OpenMP
    for (int y = 0; y < img_height; y++)
    { // Loop over image rows
        // *** Commented out for Visual Studio, fprintf is not thread-safe
        //fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
        for (unsigned short x = 0; x < img_width; x++)                               // Loop cols
            for (int sy = 0, i = (img_height - y - 1) * img_width + x; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++)
                { // 2x2 subpixel cols
                    Vector ray;
                    for (int s = 0; s < sample_per_subpixel; s++)
                    {
                        double r1 = 2 * RAND(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * RAND(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vector d = film_x_axis * (((sx + .5 + dx) / 2 + x) / img_width - .5) +
                                   film_y_axis * (((sy + .5 + dy) / 2 + y) / img_height - .5) + camera.direction;
                        ray = ray + PathTracing(Ray(camera.origin + d * 140, d.Normal()), 0) * (1. / sample_per_subpixel);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    img_buffer[i] = img_buffer[i] + Vector(ClampValue(ray.x), ClampValue(ray.y), ClampValue(ray.z)) * .25;
                }

        timer_progress.End();
        if (timer_progress.GetTime() > 1.0)
        {
            std::cout << "Rendering Progress: " << std::setiosflags(std::ios::fixed)
                      << std::setprecision(2) << (double)((y + 1) * 100.0 / img_height) << "%" << std::endl;
            timer_progress.Start();
        }
    }

    timer_rendertime.End();
    std::cout << "============" << std::endl;
    std::cout << "Finish !" << std::endl;
    std::cout << "Rednering Time: " << timer_rendertime.GetTime() << " sec" << std::endl;

    FILE *f = fopen(argv[1], "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", img_width, img_height, 255);
    for (int i = 0; i < img_width * img_height; i++)
        fprintf(f, "%d\n%d\n%d\n", ColorFloat2Int(img_buffer[i].x), ColorFloat2Int(img_buffer[i].y), ColorFloat2Int(img_buffer[i].z));
}
