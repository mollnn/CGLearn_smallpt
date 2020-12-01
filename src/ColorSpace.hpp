#ifndef _COLORSPACE_HPP
#define _COLORSPACE_HPP

#include <bits/stdc++.h>

const double pi = 3.1415926535897932384626433832795;

inline double RadianceToBrightness(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int ColorFloat2Int(double x) { return int(pow(RadianceToBrightness(x), 1 / 2.2) * 255 + .5); }

#endif