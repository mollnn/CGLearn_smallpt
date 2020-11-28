#ifndef _VECTOR_HPP
#define _VECTOR_HPP

#include <bits/stdc++.h>

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
    Vector Cross(Vector &b) {return *this%b;}
};

#endif