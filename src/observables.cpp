#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include "observables.hpp"
#include "global_decl.hpp"
#include "generic_func.hpp"
#include "lattice_ops.hpp"
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