#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>

#include "observables.hpp"
#include "simulation.hpp"
#include "lattice_ops.hpp"
#include "statistics.hpp"
double plaquette(std::unique_ptr<Simulation> &sim, int site_index, int (&shift)[4], int mu, int nu)
{
    int y[4], z[4];
    copyCoordinates(shift, z);
    copyCoordinates(shift, y);
    y[mu]++;
    z[nu]++;
    matrix temp = callLatticeSite(sim, site_index, shift, mu);
    temp *= callLatticeSite(sim, site_index, y, nu);
    temp *= callLatticeSite(sim, site_index, z, mu).adjoint();
    temp *= callLatticeSite(sim, site_index, shift, nu).adjoint();
    return 0.5 * temp.trace().real();
}

double plaquetteAverage(std::unique_ptr<Simulation> &sim)
{
    double accumulator = 0.0;
    int jumpNone[4] = {0};
    for (int site_index = 0; site_index < sim->lsites; site_index++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            for (int mu = 0; mu < nu; mu++)
            {
                accumulator += plaquette(sim, site_index, jumpNone, mu, nu);
            }
        }
    }
    return accumulator / (sim->lsites * 2 * 3);
}

double rectangle(std::unique_ptr<Simulation> &sim, int site_index, int (&shift)[4], int mu, int nu, int mu_len, int nu_len)
{
    int y[4], z[4];
    copyCoordinates(shift, z);
    copyCoordinates(shift, y);
    matrix forward = matrix::Identity(), backward = matrix::Identity();
    // Bottom
    for (int i = 0; i < mu_len; i++)
    {
        forward = forward * callLatticeSite(sim, site_index, y, mu);
        y[mu]++;
    }
    // Right
    for (int i = 0; i < nu_len; i++)
    {
        forward = forward * callLatticeSite(sim, site_index, y, nu);
        y[nu]++;
    }
    // Left
    for (int i = 0; i < nu_len; i++)
    {
        backward = backward * callLatticeSite(sim, site_index, z, nu);
        z[nu]++;
    }
    // Top
    for (int i = 0; i < mu_len; i++)
    {
        backward = backward * callLatticeSite(sim, site_index, z, mu);
        z[mu]++;
    }

    return 0.5 * (forward * backward.adjoint()).trace().real();
}

double rectangleAverage(std::unique_ptr<Simulation> &sim, int mu_len, int nu_len)
{
    double accumulator = 0.0;
    int jumpNone[4] = {0};
    for (int site_index = 0; site_index < sim->lsites; site_index++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            for (int mu = 0; mu < nu; mu++)
            {
                accumulator += rectangle(sim, site_index, jumpNone, mu, nu, mu_len, nu_len);
            }
        }
    }
    return accumulator / (sim->lsites * 2 * 3);
}

double polyakovLine(std::unique_ptr<Simulation> &sim, int site_index)
{
    int x[4];
    siteIndexToCoordinates(sim, site_index, x[0], x[1], x[2], x[3]);
    matrix temp = matrix::Identity();
    for (int i = 0; i < sim->lt; i++)
    {
        x[0] = (x[0] + 1 + sim->lt) % sim->lt;
        temp = temp * sim->lattice[coordinatesToSiteIndex(sim, x[0], x[1], x[2], x[3])].field[0];
    }
    return 0.5 * temp.trace().real();
}
double polyakovLine(std::unique_ptr<Simulation> &sim, int site_index, int dir)
{
    int x[4];
    siteIndexToCoordinates(sim, site_index, x[0], x[1], x[2], x[3]);
    matrix temp = matrix::Identity();
    for (int i = 0; i < sim->ldir[dir]; i++)
    {
        x[dir] = (x[dir] + 1 + sim->ldir[dir]) % sim->ldir[dir];
        temp = temp * sim->lattice[coordinatesToSiteIndex(sim, x[0], x[1], x[2], x[3])].field[dir];
    }
    return 0.5 * temp.trace().real();
}
void polyakovLines(std::unique_ptr<Simulation> &sim, std::string filename, int dir1, int dir2)
{
    int x[4] = {0};
    std::ofstream file(sim->datFolder + filename, std::ios_base::app);
    for (int i = 0; i < sim->ldir[dir1]; i++)
    {
        for (int j = 0; j < sim->ldir[dir2]; j++)
        {
            x[dir1] = i;
            x[dir2] = j;
            file << polyakovLine(sim, coordinatesToSiteIndex(sim, x[0], x[1], x[2], x[3])) << "\n";
        }
    }
    file.close();
}
void polyakovLinesAbs(std::unique_ptr<Simulation> &sim, std::string filename)
{
    int x[4] = {0};
    std::ofstream file(sim->datFolder + filename, std::ios_base::app);
    std::vector<double> abs;
    abs.reserve(4 * sim->l * sim->l * sim->lt);
    for (int dir = 0; dir < 4; dir++)
    {
        int dir1, dir2, dir3;
        dir1 = (dir + 1) % 4;
        dir2 = (dir + 2) % 4;
        dir3 = (dir + 3) % 4;
        x[dir] = 0;
        for (int i = 0; i < sim->ldir[dir1]; i++)
        {
            for (int j = 0; j < sim->ldir[dir2]; j++)
            {
                for (int k = 0; k < sim->ldir[dir3]; k++)
                {
                    x[dir1] = i;
                    x[dir2] = j;
                    x[dir3] = k;
                    abs.push_back(std::abs(polyakovLine(sim, coordinatesToSiteIndex(sim, x[0], x[1], x[2], x[3]), dir)));
                }
            }
        }
    }
    if (abs.size() != static_cast<long unsigned int>(4 * sim->l * sim->l * sim->lt))
        std::cout << "Wrong number of lines.\n";
    file << std::accumulate(abs.begin(), abs.end(), 0.0) / static_cast<double>(abs.size()) << std::endl;
    file.close();
}

double higgsSquareAverage(std::unique_ptr<Simulation> &sim)
{
    double accumulator = 0.0;
    for (int site_index = 0; site_index < sim->lsites; site_index++)
    {
        accumulator += (sim->lattice[site_index].field[4] * sim->lattice[site_index].field[4]).trace().real();
    }
    return accumulator / sim->lsites;
}
double higgsSquareSum(std::unique_ptr<Simulation> &sim)
{
    double accumulator = 0.0;
    matrix tmp;
    int jumpNone[4] = {0};
    for (int site_index = 0; site_index < sim->lsites; site_index++)
    {
        callLatticeSite(sim, tmp, site_index, jumpNone, 4);
        accumulator += (tmp * tmp).trace().real();
    }
    return accumulator;
}
matrix higgsAverage(std::unique_ptr<Simulation> &sim)
{
    matrix accumulator = matrix::Zero();
    for (int site_index = 0; site_index < sim->lsites; site_index++)
    {
        accumulator += (sim->lattice[site_index].field[4]);
    }
    return accumulator / sim->lsites;
}

double correlator(std::unique_ptr<Simulation> &sim, int site_index, int time_forward)
{
    matrix accumulator = sim->lattice[site_index].field[4];

    int x[4];
    siteIndexToCoordinates(sim, site_index, x[0], x[1], x[2], x[3]);

    for (int i = x[0]; i < x[0] + time_forward; i++)
    {
        accumulator = accumulator * sim->lattice[coordinatesToSiteIndex(sim, (x[0] + i + sim->ldir[0]) % sim->ldir[0], x[1], x[2], x[3])].field[0];
    }
    accumulator = accumulator * sim->lattice[coordinatesToSiteIndex(sim, (x[0] + time_forward + sim->ldir[0]) % sim->ldir[0], x[1], x[2], x[3])].field[4];

    for (int i = time_forward; i != x[0]; i = (i + 1 + sim->ldir[0]) % sim->ldir[0])
    {
        accumulator = accumulator * sim->lattice[coordinatesToSiteIndex(sim, i, x[1], x[2], x[3])].field[0];
    }
    return accumulator.trace().real();
}

double averageCorrelatorVolume(std::unique_ptr<Simulation> &sim, int time_forward)
{
    std::vector<double> data;
    for (int i = 0; i < sim->ldir[1]; i++)
    {
        for (int j = 0; j < sim->ldir[2]; j++)
        {
            for (int k = 0; k < sim->ldir[3]; k++)
            {
                data.push_back(correlator(sim, coordinatesToSiteIndex(sim, 0, i, j, k), time_forward));
            }
        }
    }
    return average(data);
}

void copyCoordinates(int (&src)[4], int (&dest)[4])
{
    dest[0] = src[0];
    dest[1] = src[1];
    dest[2] = src[2];
    dest[3] = src[3];
}
void addCoordinatesInPlace(int (&src)[4], int (&dest)[4])
{
    dest[0] += src[0];
    dest[1] += src[1];
    dest[2] += src[2];
    dest[3] += src[3];
}