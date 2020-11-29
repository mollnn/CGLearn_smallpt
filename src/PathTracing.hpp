#ifndef _PATHTRACING_HPP
#define _PATHTRACING_HPP

#include <bits/stdc++.h>

#include "FloatUtil.hpp"

#include "Timer.hpp"
#include "Random.hpp"
#include "Vector.hpp"
#include "Ray.hpp"
#include "Material.hpp"
#include "Sphere.hpp"

#include "Scene.hpp"

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

Vector PathTracing(const Ray &ray, int depth)
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
        return obj.emission + objColor.DirectMult(PathTracing(Ray(x, d), depth));
    }
    else if (obj.material == SPECULAR) // Ideal SPECULAR reflection
    {
        return obj.emission + objColor.DirectMult(PathTracing(Ray(x, ray.direction - normalOut * 2 * normalOut.Dot(ray.direction)), depth));
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
            return obj.emission + objColor.DirectMult(PathTracing(reflRay, depth));
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
            return obj.emission + objColor.DirectMult(depth > 2 ? (RAND() < P ? PathTracing(reflRay, depth) * RP
                                                                                        : PathTracing(Ray(x, tdir), depth) * TP)
                                                                : PathTracing(reflRay, depth) * Re + PathTracing(Ray(x, tdir), depth) * Tr);
        }
    }
}

#endif