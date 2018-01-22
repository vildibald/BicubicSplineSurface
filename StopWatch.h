#pragma once
#include <windows.h>


class StopWatch
{
    LARGE_INTEGER frequency_;        // ticks per second
    LARGE_INTEGER t1_, t2_;           // ticks
public:
    StopWatch()
    {
        // get ticks per second
        QueryPerformanceFrequency(&frequency_);
    }

    void Start()
    {
        QueryPerformanceCounter(&t1_);
    }

    void Stop()
    {
        QueryPerformanceCounter(&t2_);
    }

    double EllapsedTime()
    {
        return (t2_.QuadPart - t1_.QuadPart) * 1000.0
               / frequency_.QuadPart;
    }
};