#include <bits/stdc++.h>
using namespace std;


int STDRANDOM()
{
    return rand()*rand();
}

unsigned int xorshf96(void)
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

signed main()
{
    unsigned int mx=0;
    for(int i=0;i<10000;i++)
    {
        mx=max(mx,xorshf96());
    }
    cout<<mx<<endl;

}