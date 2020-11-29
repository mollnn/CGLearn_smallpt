#ifndef _SCENE_HPP
#define _SCENE_HPP

#include <bits/stdc++.h>

#include "Sphere.hpp"

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

#endif