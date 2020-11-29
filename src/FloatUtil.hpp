#ifndef _FLOATUTIL_HPP
#define _FLOATUTIL_HPP

#include <bits/stdc++.h>

const double pi = 3.1415926535897932384626433832795;

inline double ClampValue(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int ColorFloat2Int(double x) { return int(pow(ClampValue(x), 1 / 2.2) * 255 + .5); }

#endif