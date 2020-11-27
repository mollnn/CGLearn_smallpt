#include <bits/stdc++.h>
#include "Timer.hpp"
#include "Random.hpp"

#define pi 3.1415926535897932384626433832795

// 3 维向量，用于空间向量和颜色的表示
struct Vector
{
    double x, y, z;
    Vector(double x_ = 0, double y_ = 0, double z_ = 0)
    {
        x = x_;
        y = y_;
        z = z_;
    }
    Vector operator+(const Vector &b) const { return Vector(x + b.x, y + b.y, z + b.z); }
    Vector operator-(const Vector &b) const { return Vector(x - b.x, y - b.y, z - b.z); }
    Vector operator*(double b) const { return Vector(x * b, y * b, z * b); }
    Vector DirectMult(const Vector &b) const { return Vector(x * b.x, y * b.y, z * b.z); }
    Vector &Normal() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double Dot(const Vector &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
    Vector operator%(Vector &b) { return Vector(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};
struct Ray
{
    Vector origin, direction;
    Ray(Vector o_, Vector d_) : origin(o_), direction(d_) {}
};
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
Vector Radiance(const Ray &ray, int depth, unsigned short *nRandSeeds)
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
        if (RAND(nRandSeeds) < p)
            objColor = objColor * (1 / p);
        else
            return obj.emission; //R.R.
    if (obj.material == DIFFUSE) // Ideal DIFFUSE reflection
    {
        double r1 = 2 * pi * RAND(nRandSeeds);
        double r2 = RAND(nRandSeeds);
        double r2sqrt = sqrt(r2);
        Vector w = normal;
        Vector u = ((fabs(w.x) > .1 ? Vector(0, 1) : Vector(1)) % w).Normal();
        Vector v = w % u;
        Vector d = (u * cos(r1) * r2sqrt + v * sin(r1) * r2sqrt + w * sqrt(1 - r2)).Normal();
        return obj.emission + objColor.DirectMult(Radiance(Ray(x, d), depth, nRandSeeds));
    }
    else if (obj.material == SPECULAR) // Ideal SPECULAR reflection
    {
        return obj.emission + objColor.DirectMult(Radiance(Ray(x, ray.direction - normalOut * 2 * normalOut.Dot(ray.direction)), depth, nRandSeeds));
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
            return obj.emission + objColor.DirectMult(Radiance(reflRay, depth, nRandSeeds));
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
            return obj.emission + objColor.DirectMult(depth > 2 ? (RAND(nRandSeeds) < P ? Radiance(reflRay, depth, nRandSeeds) * RP
                                                                                        : Radiance(Ray(x, tdir), depth, nRandSeeds) * TP)
                                                                : Radiance(reflRay, depth, nRandSeeds) * Re + Radiance(Ray(x, tdir), depth, nRandSeeds) * Tr);
        }
    }
}
int main(int argc, char *argv[])
{
    srand(time(NULL));

    Timer timer;
    int nImageWidth = 256;
    int nImageHeight = 192;
    int nSamplePerSubpixel = 8;                                           // # samples
    Ray camera(Vector(50, 52, 295.6), Vector(0, -0.042612, -1).Normal()); // cam pos, dir
    Vector vecFilmX = Vector(nImageWidth * .5135 / nImageHeight);
    Vector vecFilmY = (vecFilmX % camera.direction).Normal() * .5135;
    Vector *imageBuffer = new Vector[nImageWidth * nImageHeight];
#pragma omp parallel for schedule(dynamic, 1) private(ray) // OpenMP
    for (int y = 0; y < nImageHeight; y++)
    { // Loop over image rows
        // *** Commented out for Visual Studio, fprintf is not thread-safe
        //fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
        unsigned short nRandSeeds[3] = {0, 0, y * y * y};                                // *** Moved outside for VS2012
        for (unsigned short x = 0; x < nImageWidth; x++)                                 // Loop cols
            for (int sy = 0, i = (nImageHeight - y - 1) * nImageWidth + x; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++)
                { // 2x2 subpixel cols
                    Vector ray;
                    for (int s = 0; s < nSamplePerSubpixel; s++)
                    {
                        double r1 = 2 * RAND(nRandSeeds), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * RAND(nRandSeeds), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vector d = vecFilmX * (((sx + .5 + dx) / 2 + x) / nImageWidth - .5) +
                                   vecFilmY * (((sy + .5 + dy) / 2 + y) / nImageHeight - .5) + camera.direction;
                        ray = ray + Radiance(Ray(camera.origin + d * 140, d.Normal()), 0, nRandSeeds) * (1. / nSamplePerSubpixel);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    imageBuffer[i] = imageBuffer[i] + Vector(ClampValue(ray.x), ClampValue(ray.y), ClampValue(ray.z)) * .25;
                }
        std::cout << "Rendering Progress: " << std::setiosflags(std::ios::fixed)
                  << std::setprecision(2) << (double)((y + 1) * 100.0 / nImageHeight) << "%" << std::endl;
    }
    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", nImageWidth, nImageHeight, 255);
    for (int i = 0; i < nImageWidth * nImageHeight; i++)
        fprintf(f, "%d\n%d\n%d\n", ColorFloat2Int(imageBuffer[i].x), ColorFloat2Int(imageBuffer[i].y), ColorFloat2Int(imageBuffer[i].z));
}
