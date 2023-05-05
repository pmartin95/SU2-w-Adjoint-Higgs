#include "statistics.hpp"
#include <numeric>
double average(const std::vector<double> &input)
{
    return std::accumulate(input.begin(), input.end(), 0.0) / input.size();
}

matrix average(const std::vector<matrix> &input)
{
    matrix sum = matrix::Zero();
    for (const auto &m : input)
    {
        sum += m;
    }
    return sum / input.size();
}

double sumVec(const std::vector<double> &input)
{
    return std::accumulate(input.begin(), input.end(), 0.0);
}

int computeJackknifeStatistics(const std::vector<double> &inputData, int setLength, double &Jackknife_ave, double &Jackknife_error)
{
    if (inputData.size() % setLength != 0)
        return JACKKNIFE_NO_DIVIDE;

    int numSets = inputData.size() / static_cast<long unsigned int>(setLength);
    std::vector<double> holdingVector;
    std::vector<double> setOfAverages(numSets, 0.0);
    for (int i = 0; i < numSets; i++)
    {
        holdingVector.clear();
        holdingVector.reserve(inputData.size() - setLength);
        holdingVector.insert(holdingVector.end(), inputData.begin(), inputData.begin() + i * setLength);
        holdingVector.insert(holdingVector.end(), inputData.begin() + (i + 1) * setLength, inputData.end());
        setOfAverages[i] = average(holdingVector);
    }

    Jackknife_ave = average(setOfAverages);

    std::vector<double> deviations(numSets, 0.0);
    double temp;
    for (int i = 0; i < numSets; i++)
    {
        temp = (setOfAverages[i] - Jackknife_ave) / Jackknife_ave;
        deviations[i] = temp * temp;
    }
    Jackknife_error = Jackknife_ave * sqrt(static_cast<double>(numSets - 1) * average(deviations));
    return JACKKNIFE_SUCCESS;
}