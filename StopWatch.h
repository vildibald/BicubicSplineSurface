#pragma once
//#include <windows.h>


class StopWatch
{
public:
    StopWatch() = default;
    std::chrono::high_resolution_clock::time_point start_;
    std::chrono::high_resolution_clock::time_point stop_;

    void Start()
    {
        start_ = std::chrono::high_resolution_clock::now();
    }

    void Stop()
    {
        stop_ = std::chrono::high_resolution_clock::now();
    }

    long long EllapsedTime()
    {
        const auto diff = stop_ - start_;
        return std::chrono::duration_cast<std::chrono::microseconds>(diff).count();
    }
};