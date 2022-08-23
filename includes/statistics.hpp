#pragma once
#include <vector>
#include <cmath>
enum jackknife
{
    JACKKNIFE_SUCCESS,
    JACKKNIFE_NO_DIVIDE,
    JACKKNIFE_UNEQUAL_SIZED_DATA
};
double average(const std::vector<double> &input);
double sumVec(const std::vector<double> &input);
int computeJackknifeStatistics(const std::vector<double> &inputData, int setLength, double &Jackknife_ave, double &Jackknife_error);
