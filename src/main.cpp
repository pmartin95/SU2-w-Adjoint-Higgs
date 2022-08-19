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
typedef Eigen::Matrix<std::complex<double>, 2, 2> matrix;
enum jackknife
{
    JACKKNIFE_SUCCESS,
    JACKKNIFE_NO_DIVIDE,
    JACKKNIFE_UNEQUAL_SIZED_DATA
};
double beta = 2.3;
int MAX_ITER = 20, iter_count;
int Naccept = 0, Nreject = 0;
int l = 8, lt = 8, lsites = l * l * l * lt;
int ldir[4] = {lt, l, l, l};
std::complex<double> I(0.0, 1.0);
double rot_size = 0.5;
typedef struct link
{
    matrix field[4];
} link;
link *lattice;
std::mt19937 rng;
std::uniform_real_distribution<double> gen(0.0, 1.0);
std::normal_distribution<double> gen_normal(0.0, 1.0);
// coordinates to site index
int coordinatesToSiteIndex(int t, int x, int y, int z);
// site index to coordinates
void siteIndexToCoordinates(int site_index, int &t, int &x, int &y, int &z);
// plaquettes
double plaquette(int site_index, int mu, int nu);
double plaquetteAverage();
double polyakovLine(int site_index);
double polyakovLine(int site_index, int dir);
void polyakovLines(std::string filename, int dir1, int dir2);
void polyakovLinesAbs(std::string filename);
double rectangle(int site_index, int mu, int nu, int mu_len, int nu_len);
double rectangleAverage(int mu_len, int nu_len);
double action();
double actionPartial(int site_index, int mu);
void generateRandomSU2(matrix &m);
void generateRandomSU2Rot(matrix &m);

void hotLattice();
void coldLattice();
double average(const std::vector<double> &input);
double sumVec(const std::vector<double> &input);
int computeJackknifeStatistics(const std::vector<double> &inputData, int setLength, double &Jackknife_ave, double &Jackknife_error);

int main(int argc, char **argv)
{
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
    beta_distribution.push_back(3.0);

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

double plaquette(int site_index, int mu, int nu)
{
    int x[4];
    siteIndexToCoordinates(site_index, x[0], x[1], x[2], x[3]);
    int y[4] = {x[0], x[1], x[2], x[3]}, z[4] = {x[0], x[1], x[2], x[3]};
    y[mu] = (x[mu] + 1) % ldir[mu];
    z[nu] = (x[nu] + 1) % ldir[nu];
    matrix temp = lattice[site_index].field[mu];

    temp *= lattice[coordinatesToSiteIndex(y[0], y[1], y[2], y[3])].field[nu];
    temp *= lattice[coordinatesToSiteIndex(z[0], z[1], z[2], z[3])].field[mu].adjoint();
    temp *= lattice[site_index].field[nu].adjoint();
    return 0.5 * temp.trace().real();
}
double plaquetteAverage()
{
    double accumulator = 0.0;
    int count = 0;
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            for (int mu = 0; mu < nu; mu++)
            {
                accumulator += plaquette(site_index, mu, nu);
            }
        }
    }
    return accumulator / (lsites * 2 * 3);
}
double rectangle(int site_index, int mu, int nu, int mu_len, int nu_len)
{
    int x[4];
    siteIndexToCoordinates(site_index, x[0], x[1], x[2], x[3]);
    int y[4] = {x[0], x[1], x[2], x[3]}, z[4] = {x[0], x[1], x[2], x[3]};

    matrix forward = matrix::Identity(), backward = matrix::Identity();

    // bottom
    for (int i = 0; i < mu_len; i++)
    {
        y[mu] = (x[mu] + i + ldir[mu]) % ldir[mu];
        forward *= lattice[coordinatesToSiteIndex(y[0], y[1], y[2], y[3])].field[mu];
    }
    y[mu] = (x[mu] + mu_len + ldir[mu]) % ldir[mu];
    // right
    for (int i = 0; i < nu_len; i++)
    {
        y[nu] = (x[nu] + i + ldir[nu]) % ldir[nu];
        forward *= lattice[coordinatesToSiteIndex(y[0], y[1], y[2], y[3])].field[nu];
    }
    // left
    for (int i = 0; i < nu_len; i++)
    {
        z[nu] = (x[nu] + i + ldir[nu]) % ldir[nu];
        backward *= lattice[coordinatesToSiteIndex(z[0], z[1], z[2], z[3])].field[nu];
    }
    z[nu] = (x[nu] + nu_len + ldir[nu]) % ldir[nu];
    // top
    for (int i = 0; i < mu_len; i++)
    {
        z[mu] = (x[mu] + i + ldir[mu]) % ldir[mu];
        backward *= lattice[coordinatesToSiteIndex(z[0], z[1], z[2], z[3])].field[mu];
    }

    return 0.5 * (forward * backward.adjoint()).trace().real();
}

double rectangleAverage(int mu_len, int nu_len)
{
    double accumulator = 0.0;
    int count = 0;
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            for (int mu = 0; mu < nu; mu++)
            {
                accumulator += rectangle(site_index, mu, nu, mu_len, nu_len);
            }
        }
    }
    return accumulator / (lsites * 2 * 3);
}

double polyakovLine(int site_index)
{
    int x[4];
    siteIndexToCoordinates(site_index, x[0], x[1], x[2], x[3]);
    matrix temp = matrix::Identity();
    for (int i = 0; i < lt; i++)
    {
        x[0] = (x[0] + 1 + lt) % lt;
        temp = temp * lattice[coordinatesToSiteIndex(x[0], x[1], x[2], x[3])].field[0];
    }
    return 0.5 * temp.trace().real();
}
double polyakovLine(int site_index, int dir)
{
    int x[4];
    siteIndexToCoordinates(site_index, x[0], x[1], x[2], x[3]);
    matrix temp = matrix::Identity();
    for (int i = 0; i < ldir[dir]; i++)
    {
        x[dir] = (x[dir] + 1 + ldir[dir]) % ldir[dir];
        temp = temp * lattice[coordinatesToSiteIndex(x[0], x[1], x[2], x[3])].field[dir];
    }
    return 0.5 * temp.trace().real();
}
void polyakovLines(std::string filename, int dir1, int dir2)
{
    int x[4] = {0};
    std::ofstream file(filename, std::ios_base::app);
    for (int i = 0; i < ldir[dir1]; i++)
    {
        for (int j = 0; j < ldir[dir2]; j++)
        {
            x[dir1] = i;
            x[dir2] = j;
            file << polyakovLine(coordinatesToSiteIndex(x[0], x[1], x[2], x[3])) << "\n";
        }
    }
    file.close();
}
void polyakovLinesAbs(std::string filename)
{
    int x[4] = {0};
    std::ofstream file(filename, std::ios_base::app);
    std::vector<double> abs;
    abs.reserve(4 * l * l * lt);
    for (int dir = 0; dir < 4; dir++)
    {
        int dir1, dir2, dir3;
        dir1 = (dir + 1) % 4;
        dir2 = (dir + 2) % 4;
        dir3 = (dir + 3) % 4;
        x[dir] = 0;
        for (int i = 0; i < ldir[dir1]; i++)
        {
            for (int j = 0; j < ldir[dir2]; j++)
            {
                for (int k = 0; k < ldir[dir3]; k++)
                {
                    x[dir1] = i;
                    x[dir2] = j;
                    x[dir3] = k;
                    abs.push_back(std::abs(polyakovLine(coordinatesToSiteIndex(x[0], x[1], x[2], x[3]), dir)));
                }
            }
        }
    }
    if (abs.size() != 4 * l * l * lt)
        std::cout << "Wrong number of lines.\n";
    file << std::accumulate(abs.begin(), abs.end(), 0.0) / abs.size() << std::endl;
    file.close();
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