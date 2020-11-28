#ifndef _RAY_HPP
#define _RAY_HPP

#include <bits/stdc++.h>
#include "Vector.hpp"

struct Ray
{
    Vector origin, direction;
    Ray(Vector o_, Vector d_) : origin(o_), direction(d_) {}
};

#endif 