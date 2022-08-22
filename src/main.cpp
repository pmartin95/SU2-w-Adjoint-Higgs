#include <complex>
#include <random>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <limits>
#include <cstdlib>

#include "global_decl.hpp"
#include "observables.hpp"
#include "generic_func.hpp"


enum jackknife
{
    JACKKNIFE_SUCCESS,
    JACKKNIFE_NO_DIVIDE,
    JACKKNIFE_UNEQUAL_SIZED_DATA
};
link *lattice;
double beta = 2.3;
int MAX_ITER = 1;
int iter_count = 0;
int Naccept = 0, Nreject = 0;
int l = 8, lt = 8, lsites = l * l * l * lt;
int ldir[4] = {lt, l, l, l};
std::complex<double> I(0.0, 1.0);
double rot_size = 0.4;

std::mt19937 rng;
std::uniform_real_distribution<double> gen(0.0, 1.0);
std::normal_distribution<double> gen_normal(0.0, 1.0);

int main(int argc, char **argv)
{
    double input_beta = 2.3;
    if (argc > 0)
    {
        input_beta = atof(argv[1]);
        std::cout << "Running beta= " << input_beta << std::endl;
    }
    // Boilerplate
    std::random_device r;
    rng.seed(r());

    lattice = new link[lsites];
    matrix tmp;
    double actionDiff;

    // File setup
    int highest_len = 3;
    std::string base_filename = "rect";
    std::vector<std::string> MxMfilenames, MxNfilenames;
    std::vector<std::ofstream> MxMfiles, MxNfiles;
    std::vector<double> rect2x2observe;
    for (int i = 0; i < highest_len; i++) // Name files and open them
    {
        MxMfilenames.push_back(base_filename + std::to_string(i + 1) + "x" + std::to_string(i + 1) + ".txt");
        MxNfilenames.push_back(base_filename + std::to_string(i + 1) + "x" + std::to_string(i) + ".txt");
        std::ofstream temp1, temp2;
        temp1.open(MxMfilenames[i]);
        temp2.open(MxNfilenames[i]);
        MxMfiles.push_back(std::move(temp1));
        MxNfiles.push_back(std::move(temp2));
    }

    // Setup Beta distribution 0.0 to 4.0 with 0.1 in between
    std::vector<double> beta_distribution; // Create beta distribution
    // {
    //     double temp = 0.1;
    //     do
    //     {
    //         beta_distribution.push_back(temp);
    //         temp += 0.1;
    //     } while (temp < 4.0);
    // }
    beta_distribution.push_back(input_beta);

    for (auto &B : beta_distribution)
    {
        beta = B;
        std::vector<std::vector<double>> MxMdata, MxNdata;
        for (int i = 0; i < highest_len; i++) // Init data holding vectors
        {
            std::vector<double> tempdata1, tempdata2;
            MxMdata.push_back(std::move(tempdata1));
            MxNdata.push_back(std::move(tempdata2));
        }
        std::vector<double> MxMave(highest_len, 0.0), MxMerror(highest_len, 0.0), MxNave(highest_len, 0.0), MxNerror(highest_len, 0.0);
        double ave, error;
        // Thermalization
        hotLattice();
        for (int i = 0; i < 5000; i++)
        {
            for (int site_index = 0; site_index < lsites; site_index++)
            {

                for (int dir = 0; dir < 4; dir++)
                {

                    actionDiff = actionPartial(site_index, dir);
                    tmp = lattice[site_index].field[dir];
                    generateRandomSU2Rot(lattice[site_index].field[dir]);
                    actionDiff -= actionPartial(site_index, dir);
                    // actionDiff = -Delta S
                    actionDiff = std::exp(actionDiff);
                    if (actionDiff < gen(rng))
                    {
                        lattice[site_index].field[dir] = tmp;
                        ++Nreject;
                    }
                    else
                    {
                        ++Naccept;
                    }
                }
            }
            polyakovLinesAbs("hotAbsPoly" + std::to_string(beta) + ".txt");
            // std::cout << rectangleAverage(2, 1) << std::endl;
            // std::cout << lattice[5].field[2] << std::endl;
            std::cout << "Acceptance rate in thermalization:" << static_cast<double>(Naccept) / static_cast<double>(Naccept + Nreject) << std::endl;
            Naccept = 0;
            Nreject = 0;
        }

        iter_count = 0;
        do
        {
            ++iter_count;
            for (int i = 0; i < 100; i++)
            {
                for (int j = 0; j < 1000; j++)
                {
                    for (int site_index = 0; site_index < lsites; site_index++)
                    {

                        for (int dir = 0; dir < 4; dir++)
                        {
                            actionDiff = actionPartial(site_index, dir);
                            tmp = lattice[site_index].field[dir];
                            generateRandomSU2Rot(lattice[site_index].field[dir]);
                            actionDiff -= actionPartial(site_index, dir); // actionDiff = -Delta S

                            actionDiff = std::exp(actionDiff);
                            if (actionDiff < gen(rng))
                            {
                                lattice[site_index].field[dir] = tmp;
                                ++Nreject;
                            }
                            else
                            {
                                ++Naccept;
                            }
                        }
                    }
                    // std::cout << "Acceptance: " << static_cast<double>(Naccept) / static_cast<double>(Naccept + Nreject) << std::endl;
                    Naccept = 0;
                    Nreject = 0;
                }
                for (int i = 0; i < highest_len; i++)
                {
                    MxMdata[i].push_back(rectangleAverage(i + 1, i + 1));
                    MxNdata[i].push_back(rectangleAverage(i + 1, i));
                }
                polyakovLines("polyakov" + std::to_string(beta) + ".txt", 1, 2);
            }
            for (int i = 0; i < highest_len; i++)
            {
                computeJackknifeStatistics(MxMdata[i], 10, MxMave[i], MxMerror[i]);
                computeJackknifeStatistics(MxNdata[i], 10, MxNave[i], MxNerror[i]);
            }
            std::cout << "beta: " << beta << " X(1,1) = " << MxMave[0] << "\u00b1" << MxMerror[0] << " X(2,2)%=" << std::abs(MxMerror[1] / MxMave[1]) << " X(1,2)%=" << std::abs(MxNerror[1] / MxNave[1]) << std::endl;
        } while ((std::abs(MxMerror[0] / MxMave[0]) > 0.05 || std::abs(MxMerror[1] / MxMave[1]) > 0.1 || std::abs(MxNerror[1] / MxNave[1]) > 0.1) && (iter_count < MAX_ITER));

        for (int i = 0; i < highest_len; i++)
        {
            MxMfiles[i] << beta << " " << MxMave[i] << " " << MxMerror[i] << "\n";
            MxNfiles[i] << beta << " " << MxNave[i] << " " << MxNerror[i] << "\n";
        }
    }

    delete[] lattice;
    for (auto &file : MxMfiles)
        file.close();
    for (auto &file : MxNfiles)
        file.close();

    return 0;
}

int coordinatesToSiteIndex(int t, int x, int y, int z)
{
    int site_index = t;
    site_index = site_index * l + x;
    site_index = site_index * l + y;
    site_index = site_index * l + z;

    return site_index;
}

void siteIndexToCoordinates(int site_index, int &t, int &x, int &y, int &z)
{
    int temp = site_index;
    z = temp % l;
    temp = temp / l;
    y = temp % l;
    temp = temp / l;
    x = temp % l;
    t = temp / l;
}
double action()
{
    double accumulator = 0.0;
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            for (int mu = 0; mu < nu; mu++)
            {
                accumulator += (1.0 - plaquette(site_index, mu, nu));
            }
        }
    }
    return beta * accumulator;
}
double actionPartial(int site_index, int mu)
{
    double accumulator = 0.0;
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu)
            continue;
        int x[4];
        siteIndexToCoordinates(site_index, x[0], x[1], x[2], x[3]);
        x[nu] = (x[nu] - 1 + ldir[nu]) % ldir[nu];
        accumulator += plaquette(site_index, mu, nu);
        accumulator += plaquette(coordinatesToSiteIndex(x[0], x[1], x[2], x[3]), mu, nu);
    }
    return -beta * accumulator;
}

void generateRandomSU2(matrix &m)
{
    double a, b, c, d;
    double mag;
    do
    {
        a = gen_normal(rng);
        b = gen_normal(rng);
        c = gen_normal(rng);
        d = gen_normal(rng);
        mag = std::sqrt(a * a + b * b + c * c + d * d);
    } while (std::numeric_limits<double>::epsilon() * 100.0 > mag);

    m << a + b * I, c + d * I, -c + d * I, a - b * I;
    m = m / mag;
}
void generateRandomSU2Rot(matrix &m)
{
    matrix temp;
    double a, b, c, d, mag, mag_rot, rot = rot_size * gen(rng);
    do
    {
        b = gen_normal(rng);
        c = gen_normal(rng);
        d = gen_normal(rng);
        mag = std::sqrt(b * b + c * c + d * d);
    } while (std::numeric_limits<double>::epsilon() * 100.0 > mag);
    mag_rot = rot / mag;
    b = b * mag_rot;
    c = c * mag_rot;
    d = d * mag_rot;
    a = std::sqrt(1.0 - rot * rot);
    temp
        << a + b * I,
        c + d * I, -c + d * I, a - b * I;

    m = temp * m;
}
// void generateRandomSU2Rot(matrix &m)
// {
//     matrix temp;
//     generateRandomSU2(temp);
//     m = temp * m;
// }

void hotLattice()
{
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            generateRandomSU2(lattice[site_index].field[dir]);
        }
    }
}

void coldLattice()
{
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            lattice[site_index].field[dir] = matrix::Identity();
        }
    }
}

double average(const std::vector<double> &input)
{
    double temp;
    temp = sumVec(input) / static_cast<double>(input.size());
    return temp;
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
