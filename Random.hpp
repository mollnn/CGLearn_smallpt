#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#include <bits/stdc++.h>

int STDRANDOM()
{
    return rand()*rand();
}

unsigned int xorshf96()
{
    static unsigned int x = STDRANDOM(), y = STDRANDOM(), z = STDRANDOM();
    unsigned int t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;
    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;
    return z;
}

// 随机数生成器
double RAND(unsigned short seed[3])
{
    return (double)xorshf96() / (double)0xffffffffU;
}

#endif