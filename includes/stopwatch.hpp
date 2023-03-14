#include <cstdlib>
#include <sys/time.h>
#pragma once

class Timer
{
public:
    void stopwatchStart();
    double stopwatchReadSeconds();

private:
    struct timeval myStartTime;
};