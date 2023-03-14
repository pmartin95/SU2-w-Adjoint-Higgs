#include "stopwatch.hpp"
void Timer::stopwatchStart()
{
    gettimeofday(&myStartTime, NULL);
}

double Timer::stopwatchReadSeconds()
{
    struct timeval endTime;
    gettimeofday(&endTime, NULL);

    long ds = endTime.tv_sec - myStartTime.tv_sec;
    long dus = endTime.tv_usec - myStartTime.tv_usec;
    return ds + 0.000001 * dus;
}