#include <bits/stdc++.h>
#include "Timer.hpp"
#include "Random.hpp"

#define pi 3.1415926535897932384626433832795

// 3 维向量，用于空间向量和颜色的表示
struct Vec
{
    double x, y, z;
    Vec(double x_ = 0, double y_ = 0, double z_ = 0)
    {
        x = x_;
        y = y_;
        z = z_;
    }
    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec DirectMult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec &Normal() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double Dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
    Vec operator%(Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};
struct Ray
{
    Vec origin, direction;
    Ray(Vec o_, Vec d_) : origin(o_), direction(d_) {}
};
enum Refl_t
{
    DIFFUSE,
    SPECULAR,
    REFRECT
}; // material types, used in radiance()
struct Sphere
{
    double radius; // radius
    Vec p, e, c;   // position, emission, color
    Refl_t refl;   // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : radius(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &ray) const
    {                            // returns distance, 0 if nohit
        Vec op = p - ray.origin; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
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
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFFUSE),   //Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFFUSE), //Rght
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFFUSE),         //Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFFUSE),               //Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFFUSE),         //Botm
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFFUSE), //Top
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPECULAR),       //Mirr
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFRECT),        //Glas
    Sphere(600, Vec(50, 681.6 - .27, 21.6), Vec(12, 12, 12), Vec(), DIFFUSE)     //Lite
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
Vec Radiance(const Ray &ray, int depth, unsigned short *nRandSeeds)
{
    double dist; // distance to intersection
    int id = 0;  // id of intersected object
    if (!GetIntersect(ray, dist, id))
        return Vec();                                                           // if miss, return black
    const Sphere &obj = spheres[id];                                            // the hit object
    Vec x = ray.origin + ray.direction * dist;                                  // intersection point
    Vec normalOut = (x - obj.p).Normal();                                       // out-normal of sphere
    Vec normal = normalOut.Dot(ray.direction) < 0 ? normalOut : normalOut * -1; // fixed surface normal
    Vec objColor = obj.c;                                                       // object color
    double p = objColor.x > objColor.y && objColor.x > objColor.z
                   ? objColor.x
                   : objColor.y > objColor.z
                         ? objColor.y
                         : objColor.z; // max reflection
    if (depth > 100)
        return obj.e; //  prevent stack overflow
    if (++depth > 5)
        if (RAND(nRandSeeds) < p)
            objColor = objColor * (1 / p);
        else
            return obj.e;    //R.R.
    if (obj.refl == DIFFUSE) // Ideal DIFFUSE reflection
    {
        double r1 = 2 * pi * RAND(nRandSeeds);
        double r2 = RAND(nRandSeeds);
        double r2sqrt = sqrt(r2);
        Vec w = normal;
        Vec u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).Normal();
        Vec v = w % u;
        Vec d = (u * cos(r1) * r2sqrt + v * sin(r1) * r2sqrt + w * sqrt(1 - r2)).Normal();
        return obj.e + objColor.DirectMult(Radiance(Ray(x, d), depth, nRandSeeds));
    }
    else if (obj.refl == SPECULAR) // Ideal SPECULAR reflection
    {
        return obj.e + objColor.DirectMult(Radiance(Ray(x, ray.direction - normalOut * 2 * normalOut.Dot(ray.direction)), depth, nRandSeeds));
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
            return obj.e + objColor.DirectMult(Radiance(reflRay, depth, nRandSeeds));
        }
        else
        {
            Vec tdir = (ray.direction * nnt - normalOut * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).Normal();
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
            return obj.e + objColor.DirectMult(depth > 2 ? (RAND(nRandSeeds) < P ? Radiance(reflRay, depth, nRandSeeds) * RP
                                                                                 : Radiance(Ray(x, tdir), depth, nRandSeeds) * TP)
                                                         : Radiance(reflRay, depth, nRandSeeds) * Re + Radiance(Ray(x, tdir), depth, nRandSeeds) * Tr);
        }
    }
}
int main(int argc, char *argv[])
{
    srand(time(NULL));

    Timer timer;
    int w = 256, h = 192, samps = 8;                             // # samples
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).Normal()); // cam pos, dir
    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.direction).Normal() * .5135, ray, *c = new Vec[w * h];
#pragma omp parallel for schedule(dynamic, 1) private(ray) // OpenMP
    for (int y = 0; y < h; y++)
    { // Loop over image rows
        // *** Commented out for Visual Studio, fprintf is not thread-safe
        //fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
        unsigned short nRandSeeds[3] = {0, 0, y * y * y};           // *** Moved outside for VS2012
        for (unsigned short x = 0; x < w; x++)                      // Loop cols
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, ray = Vec())
                { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * RAND(nRandSeeds), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * RAND(nRandSeeds), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.direction;
                        ray = ray + Radiance(Ray(cam.origin + d * 140, d.Normal()), 0, nRandSeeds) * (1. / samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(ClampValue(ray.x), ClampValue(ray.y), ClampValue(ray.z)) * .25;
                }
        std::cout << "Rendering Progress: " << std::setiosflags(std::ios::fixed)
                  << std::setprecision(2) << (double)((y + 1) * 100.0 / h) << "%" << std::endl;
    }
    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d\n%d\n%d\n", ColorFloat2Int(c[i].x), ColorFloat2Int(c[i].y), ColorFloat2Int(c[i].z));
}
