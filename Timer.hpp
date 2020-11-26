#ifndef _TIMER_HPP
#define _TIMER_HPP

#include <bits/stdc++.h>
#include<windows.h>
#pragma comment( lib,"winmm.lib" )

struct Timer
{
    LARGE_INTEGER t1, t2, tc;
    void Start()//开始计时
    {
        QueryPerformanceFrequency(&tc);
        QueryPerformanceCounter(&t1);
    }
    void End()//结束计时
    {
        QueryPerformanceCounter(&t2);
    }
    Timer()
    {
        Start();
    }
    double GetTime()// 获取计时结果
    {
        double ans = (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart;
        return ans;
    }
    void Print()// 直接打印计时结果
    {
        double ans = (t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart;
        printf("%.6f", ans);
    }
};

#endif 