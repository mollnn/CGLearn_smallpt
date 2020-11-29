#ifndef _SPHERE_HPP
#define _SPHERE_HPP

#include <bits/stdc++.h>
#include "Vector.hpp"
#include "Ray.hpp"
#include "Material.hpp"

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

#endif