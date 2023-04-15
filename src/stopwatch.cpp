#include "stopwatch.hpp"

void Timer::stopwatchStart()
{
    myStartTime = std::chrono::high_resolution_clock::now();
}

double Timer::stopwatchReadSeconds()
{
    auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - myStartTime);
        return duration.count() * 1e-6;
}

void TimerOld::stopwatchStart()
{
    gettimeofday(&myStartTime, NULL);
}

double TimerOld::stopwatchReadSeconds()
{
    struct timeval endTime;
    gettimeofday(&endTime, NULL);

    long ds = endTime.tv_sec - myStartTime.tv_sec;
    long dus = endTime.tv_usec - myStartTime.tv_usec;
    return ds + 0.000001 * dus;
}

