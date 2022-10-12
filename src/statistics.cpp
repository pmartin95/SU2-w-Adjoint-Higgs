#include "statistics.hpp"

double average(const std::vector<double> &input)
{
    double temp;
    temp = sumVec(input) / static_cast<double>(input.size());
    return temp;
}

matrix average(const std::vector<matrix> &input)
{
    matrix temp = matrix::Zero();
    for (auto &x : input)
    {
        temp += x;
    }
    return temp / input.size();
}

double sumVec(const std::vector<double> &input)
{
    double temp_sum = 0.0;
    int vec_size = input.size();

    for (int i = 0; i < vec_size; i++)
        temp_sum += input[i];
    return temp_sum;
}

int computeJackknifeStatistics(const std::vector<double> &inputData, int setLength, double &Jackknife_ave, double &Jackknife_error)
{
    if (inputData.size() % setLength != 0)
        return JACKKNIFE_NO_DIVIDE;
    std::vector<double> holdingVector;

    int numSets;
    numSets = inputData.size() / setLength;

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