#include <bits/stdc++.h>
#include "Timer.hpp"
#include "Random.hpp"
#include "Vector.hpp"
#include "Ray.hpp"

#define pi 3.1415926535897932384626433832795

enum Material
{
    DIFFUSE,
    SPECULAR,
    REFRECT
}; // material types, used in radiance()
struct Sphere
{
    double radius;                    // radius
    Vector position, emission, color; // position, emission, color
    Material material;                // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vector p_, Vector e_, Vector c_, Material refl_) : radius(rad_), position(p_), emission(e_), color(c_), material(refl_) {}
    double intersect(const Ray &ray) const
    {                                      // returns distance, 0 if nohit
        Vector op = position - ray.origin; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps = 1e-4, b = op.Dot(ray.direction), det = b * b - op.Dot(op) + radius * radius;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};
Sphere spheres[] = {
    //Scene: radius, position, emission, color, material
    Sphere(1e5, Vector(1e5 + 1, 40.8, 81.6), Vector(), Vector(.75, .25, .25), DIFFUSE),   //Left
    Sphere(1e5, Vector(-1e5 + 99, 40.8, 81.6), Vector(), Vector(.25, .25, .75), DIFFUSE), //Rght
    Sphere(1e5, Vector(50, 40.8, 1e5), Vector(), Vector(.75, .75, .75), DIFFUSE),         //Back
    Sphere(1e5, Vector(50, 40.8, -1e5 + 170), Vector(), Vector(), DIFFUSE),               //Frnt
    Sphere(1e5, Vector(50, 1e5, 81.6), Vector(), Vector(.75, .75, .75), DIFFUSE),         //Botm
    Sphere(1e5, Vector(50, -1e5 + 81.6, 81.6), Vector(), Vector(.75, .75, .75), DIFFUSE), //Top
    Sphere(16.5, Vector(27, 16.5, 47), Vector(), Vector(1, 1, 1) * .999, SPECULAR),       //Mirr
    Sphere(16.5, Vector(73, 16.5, 78), Vector(), Vector(1, 1, 1) * .999, REFRECT),        //Glas
    Sphere(600, Vector(50, 681.6 - .27, 21.6), Vector(12, 12, 12), Vector(), DIFFUSE)     //Lite
};
inline double ClampValue(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int ColorFloat2Int(double x) { return int(pow(ClampValue(x), 1 / 2.2) * 255 + .5); }
inline bool GetIntersect(const Ray &ray, double &t, int &id)
{
    double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
    for (int i = int(n); i--;)
        if ((d = spheres[i].intersect(ray)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < inf;
}
Vector Radiance(const Ray &ray, int depth)
{
    double dist; // distance to intersection
    int id = 0;  // id of intersected object
    if (!GetIntersect(ray, dist, id))
        return Vector();                                                           // if miss, return black
    const Sphere &obj = spheres[id];                                               // the hit object
    Vector x = ray.origin + ray.direction * dist;                                  // intersection point
    Vector normalOut = (x - obj.position).Normal();                                // out-normal of sphere
    Vector normal = normalOut.Dot(ray.direction) < 0 ? normalOut : normalOut * -1; // fixed surface normal
    Vector objColor = obj.color;                                                   // object color
    double p = objColor.x > objColor.y && objColor.x > objColor.z
                   ? objColor.x
                   : objColor.y > objColor.z
                         ? objColor.y
                         : objColor.z; // max reflection
    if (depth > 100)
        return obj.emission; //  prevent stack overflow
    if (++depth > 5)
        if (RAND() < p)
            objColor = objColor * (1 / p);
        else
            return obj.emission; //R.R.
    if (obj.material == DIFFUSE) // Ideal DIFFUSE reflection
    {
        double r1 = 2 * pi * RAND();
        double r2 = RAND();
        double r2sqrt = sqrt(r2);
        Vector w = normal;
        Vector u = ((fabs(w.x) > .1 ? Vector(0, 1) : Vector(1)) % w).Normal();
        Vector v = w % u;
        Vector d = (u * cos(r1) * r2sqrt + v * sin(r1) * r2sqrt + w * sqrt(1 - r2)).Normal();
        return obj.emission + objColor.DirectMult(Radiance(Ray(x, d), depth));
    }
    else if (obj.material == SPECULAR) // Ideal SPECULAR reflection
    {
        return obj.emission + objColor.DirectMult(Radiance(Ray(x, ray.direction - normalOut * 2 * normalOut.Dot(ray.direction)), depth));
    }
    else
    {
        Ray reflRay(x, ray.direction - normalOut * 2 * normalOut.Dot(ray.direction)); // Ideal dielectric REFRACTION
        bool into = normalOut.Dot(normal) > 0;                                        // Ray from outside going in?
        double nc = 1;
        double nt = 1.5;
        double nnt = into ? nc / nt : nt / nc;
        double ddn = ray.direction.Dot(normal);
        double cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
        {
            return obj.emission + objColor.DirectMult(Radiance(reflRay, depth));
        }
        else
        {
            Vector tdir = (ray.direction * nnt - normalOut * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).Normal();
            double a = nt - nc;
            double b = nt + nc;
            double R0 = a * a / (b * b);
            double c = 1 - (into ? -ddn : tdir.Dot(normalOut));
            double Re = R0 + (1 - R0) * c * c * c * c * c;
            double Tr = 1 - Re;
            double P = .25 + .5 * Re;
            double RP = Re / P;
            double TP = Tr / (1 - P);
            // Russian roulette
            return obj.emission + objColor.DirectMult(depth > 2 ? (RAND() < P ? Radiance(reflRay, depth) * RP
                                                                                        : Radiance(Ray(x, tdir), depth) * TP)
                                                                : Radiance(reflRay, depth) * Re + Radiance(Ray(x, tdir), depth) * Tr);
        }
    }
}
int main(int argc, char *argv[])
{
    srand(time(NULL));

    Timer timer;
    int img_width = 256;
    int img_height = 192;
    int sample_per_subpixel = 8;                                           // # samples
    Ray camera(Vector(50, 52, 295.6), Vector(0, -0.042612, -1).Normal()); // cam pos, dir
    Vector film_x_axis = Vector(img_width * .5135 / img_height);
    Vector film_y_axis = (film_x_axis % camera.direction).Normal() * .5135;
    Vector *img_buffer = new Vector[img_width * img_height];
#pragma omp parallel for schedule(dynamic, 1) private(ray) // OpenMP
    for (int y = 0; y < img_height; y++)
    { // Loop over image rows
        // *** Commented out for Visual Studio, fprintf is not thread-safe
        //fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
        for (unsigned short x = 0; x < img_width; x++)                                 // Loop cols
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
                        ray = ray + Radiance(Ray(camera.origin + d * 140, d.Normal()), 0) * (1. / sample_per_subpixel);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    img_buffer[i] = img_buffer[i] + Vector(ClampValue(ray.x), ClampValue(ray.y), ClampValue(ray.z)) * .25;
                }
        std::cout << "Rendering Progress: " << std::setiosflags(std::ios::fixed)
                  << std::setprecision(2) << (double)((y + 1) * 100.0 / img_height) << "%" << std::endl;
    }
    FILE *f = fopen(argv[1], "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", img_width, img_height, 255);
    for (int i = 0; i < img_width * img_height; i++)
        fprintf(f, "%d\n%d\n%d\n", ColorFloat2Int(img_buffer[i].x), ColorFloat2Int(img_buffer[i].y), ColorFloat2Int(img_buffer[i].z));
}
