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
    if (hitobj.material == DIFFUSE) // 漫反射
    {
        // 随机生成一个光线方向
        double phi = 2 * pi * RAND();
        double radius_pow2 = RAND();
        double radius = sqrt(radius_pow2);
        // 生成微平面上的坐标轴
        Vector unitvec_z = normal;
        Vector unitvec_x = ((fabs(unitvec_z.x) > .1 ? Vector(0, 1) : Vector(1)) % unitvec_z).Normal();
        Vector unitvec_y = unitvec_z % unitvec_x;
        // 计算漫反射光线方向
        Vector diffuse_ray_dir = (unitvec_x * cos(phi) * radius + unitvec_y * sin(phi) * radius + unitvec_z * sqrt(1 - radius_pow2)).Normal();
        // 递归计算漫反射
        return hitobj.emission + hitobj_color.DirectMult(PathTracing(Ray(hit_point, diffuse_ray_dir), depth));
    }
    else if (hitobj.material == SPECULAR) // 高光反射
    {
        // 递归计算高光反射
        return hitobj.emission + hitobj_color.DirectMult(PathTracing(Ray(hit_point, ray.direction - normal_out * 2 * normal_out.Dot(ray.direction)), depth));
    }
    else if (hitobj.material == REFRECT) // 折射
    {
        Ray reflect_ray(hit_point, ray.direction - normal_out * 2 * normal_out.Dot(ray.direction));                                // 反射光线
        bool flag_into_sphere = normal_out.Dot(normal) > 0;                                                                        // 光线的方向
        double refrect_index_out = 1;                                                                                              // 外部折射率
        double refrect_index_in = 1.5;                                                                                             // 内部折射率
        double refrect_index_eff = flag_into_sphere ? refrect_index_out / refrect_index_in : refrect_index_in / refrect_index_out; // 相对折射率
        double cos_i = ray.direction.Dot(normal);                                                                                  // cos 入射角
        double cos_t_pow2 = 1 - refrect_index_eff * refrect_index_eff * (1 - cos_i * cos_i);
        if (cos_t_pow2 < 0) // 发生全反射
        {
            // 仅递归计算反射光线
            return hitobj.emission + hitobj_color.DirectMult(PathTracing(reflect_ray, depth));
        }
        else
        {
            Vector refrect_dir = (ray.direction * refrect_index_eff - normal_out * ((flag_into_sphere ? 1 : -1) * (cos_i * refrect_index_eff + sqrt(cos_t_pow2)))).Normal(); // 折射光线
            double refrect_index_delta = refrect_index_in - refrect_index_out;                                                                                               // 折射率之差
            double refrect_index_sum = refrect_index_in + refrect_index_out;                                                                                                 // 折射率之和
            double fresnel_i0 = refrect_index_delta * refrect_index_delta / (refrect_index_sum * refrect_index_sum);                                                         // 菲涅尔反射强度项
            double fresnel_tmp = 1 - (flag_into_sphere ? -cos_i : refrect_dir.Dot(normal_out));                                                                              // 用于计算实际反射强度的临时项，与角度相关
            double reflect_intensity = fresnel_i0 + (1 - fresnel_i0) * fresnel_tmp * fresnel_tmp * fresnel_tmp * fresnel_tmp * fresnel_tmp;                                  // 反射强度
            double refrect_intensity = 1 - reflect_intensity;                                                                                                                // 折射强度
            double reflect_probability = .25 + .5 * reflect_intensity;                                                                                                       // 生成反射光线的概率
            double refrect_probability = 1 - reflect_probability;                                                                                                            // 生成折射光线的概率
            double reflect_coefficient = reflect_intensity / reflect_probability;                                                                                            // 反射光线的权重系数
            double refrect_coefficient = refrect_intensity / refrect_probability;                                                                                            // 折射光线的权重系数
            // Russian roulette
            return hitobj.emission + hitobj_color.DirectMult(depth > 2 ? (RAND() < reflect_probability ? PathTracing(reflect_ray, depth) * reflect_coefficient
                                                                                                       : PathTracing(Ray(hit_point, refrect_dir), depth) * refrect_coefficient)
                                                                       : PathTracing(reflect_ray, depth) * reflect_intensity + PathTracing(Ray(hit_point, refrect_dir), depth) * refrect_intensity);
        }
    }
}

#endif