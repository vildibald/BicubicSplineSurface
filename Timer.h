//
// Created by Viliam on 28.09.2017.
//

#ifndef SPLINEKNOTS_TIMER_H
#define SPLINEKNOTS_TIMER_H

#include "StopWatch.h"

class Timer final {
    StopWatch executionWatch;
    StopWatch allWatch;
public:
    double ExecutionTime() {
        return executionWatch.EllapsedTime();
    }

    double AllTime() {
        return allWatch.EllapsedTime();
    }

    void StartExecutionTime() {
        executionWatch.Start();
    }

    void StartAllTime() {
        allWatch.Start();
    }

    void StopExecutionTime() {
        executionWatch.Stop();
    }

    void StopAllTime() {
        allWatch.Stop();
    }

    void Start() {
        allWatch.Start();
        executionWatch.Start();
    }

    void Stop() {
        allWatch.Stop();
        executionWatch.Stop();
    }

    void Reset() {
        executionWatch = StopWatch();
        allWatch = StopWatch();
    }
};
#endif //SPLINEKNOTS_TIMER_H
