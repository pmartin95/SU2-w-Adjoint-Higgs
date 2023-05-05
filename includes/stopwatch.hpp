#include <cstdlib>
#include <sys/time.h>
#include <chrono>
#pragma once

class Timer
{
public:
    void stopwatchStart();
    double stopwatchReadSeconds();

private:
    std::chrono::high_resolution_clock::time_point myStartTime;
};

class TimerOld
{
public:
    void stopwatchStart();
    double stopwatchReadSeconds();

private:
    struct timeval myStartTime;
};
