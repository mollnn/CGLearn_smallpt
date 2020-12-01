#ifndef _PATHTRACING_HPP
#define _PATHTRACING_HPP

#include <bits/stdc++.h>

#include "ColorSpace.hpp"

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
    double dist_ray_hitpoint; // 光线起点到命中点的距离
    int hitobj_id = 0;        // 命中物体编号
    // 光线与场景物体求交，若未命中则返回背景色
    if (!GetIntersect(ray, dist_ray_hitpoint, hitobj_id))
        return Vector();
    const Sphere &hitobj = spheres[hitobj_id];                                        // 命中物体
    Vector hit_point = ray.origin + ray.direction * dist_ray_hitpoint;                // 命中点
    Vector normal_out = (hit_point - hitobj.position).Normal();                       // 命中球体的外向法线
    Vector normal = normal_out.Dot(ray.direction) < 0 ? normal_out : normal_out * -1; // 实际发生折反射的法线
    Vector hitobj_color = hitobj.color;                                               // 命中物体的颜色
    double p = hitobj_color.x > hitobj_color.y && hitobj_color.x > hitobj_color.z
                   ? hitobj_color.x
                   : hitobj_color.y > hitobj_color.z
                         ? hitobj_color.y
                         : hitobj_color.z; // 最大反射强度，用于退出判断
    if (depth > 100)
        return hitobj.emission; // prevent stack overflow
    // 俄罗斯轮盘
    if (++depth > 5)
        if (RAND() < p)
            hitobj_color = hitobj_color * (1 / p);
        else
            return hitobj.emission; //R.R.
    if (hitobj.material == DIFFUSE) // Ideal DIFFUSE reflection
    {
        double phi = 2 * pi * RAND();
        double radius_pow2 = RAND();
        double radius = sqrt(radius_pow2);
        Vector unitvec_z = normal;
        Vector unitvec_x = ((fabs(unitvec_z.x) > .1 ? Vector(0, 1) : Vector(1)) % unitvec_z).Normal();
        Vector unitvec_y = unitvec_z % unitvec_x;
        Vector diffuse_ray_dir = (unitvec_x * cos(phi) * radius + unitvec_y * sin(phi) * radius + unitvec_z * sqrt(1 - radius_pow2)).Normal();
        return hitobj.emission + hitobj_color.DirectMult(PathTracing(Ray(hit_point, diffuse_ray_dir), depth));
    }
    else if (hitobj.material == SPECULAR) // Ideal SPECULAR reflection
    {
        return hitobj.emission + hitobj_color.DirectMult(PathTracing(Ray(hit_point, ray.direction - normal_out * 2 * normal_out.Dot(ray.direction)), depth));
    }
    else if (hitobj.material == REFRECT)
    {
        Ray reflect_ray(hit_point, ray.direction - normal_out * 2 * normal_out.Dot(ray.direction)); // Ideal dielectric REFRACTION
        bool flag_into_sphere = normal_out.Dot(normal) > 0;                                         // Ray from outside going in?
        double refrect_index_out = 1;
        double refrect_index_in = 1.5;
        double refrect_index_eff = flag_into_sphere ? refrect_index_out / refrect_index_in : refrect_index_in / refrect_index_out;
        double cos_i = ray.direction.Dot(normal);
        double cos_2t = 1 - refrect_index_eff * refrect_index_eff * (1 - cos_i * cos_i);
        if (cos_2t < 0) // Total internal reflection
        {
            return hitobj.emission + hitobj_color.DirectMult(PathTracing(reflect_ray, depth));
        }
        else
        {
            Vector refrect_dir = (ray.direction * refrect_index_eff - normal_out * ((flag_into_sphere ? 1 : -1) * (cos_i * refrect_index_eff + sqrt(cos_2t)))).Normal();
            double refrect_index_delta = refrect_index_in - refrect_index_out;
            double refrect_index_sum = refrect_index_in + refrect_index_out;
            double fresnel_i0 = refrect_index_delta * refrect_index_delta / (refrect_index_sum * refrect_index_sum);
            double res_cos = 1 - (flag_into_sphere ? -cos_i : refrect_dir.Dot(normal_out));
            double reflect_intensity = fresnel_i0 + (1 - fresnel_i0) * res_cos * res_cos * res_cos * res_cos * res_cos;
            double refrect_intensity = 1 - reflect_intensity;
            double reflect_probability = .25 + .5 * reflect_intensity;
            double reflect_coefficient = reflect_intensity / reflect_probability;
            double refrect_coefficient = refrect_intensity / (1 - reflect_probability);
            // Russian roulette
            return hitobj.emission + hitobj_color.DirectMult(depth > 2 ? (RAND() < reflect_probability ? PathTracing(reflect_ray, depth) * reflect_coefficient
                                                                                     : PathTracing(Ray(hit_point, refrect_dir), depth) * refrect_coefficient)
                                                                       : PathTracing(reflect_ray, depth) * reflect_intensity + PathTracing(Ray(hit_point, refrect_dir), depth) * refrect_intensity);
        }
    }
}

#endif