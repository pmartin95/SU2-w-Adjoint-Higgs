#include <complex>
#include <random>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <string.h>
#include <vector>
#include <limits>
#include <cstdlib>

#include "global_decl.hpp"
#include "observables.hpp"
#include "generic_func.hpp"
#include "metropolis-hastings.hpp"
#include "statistics.hpp"
#include "configuration_io.hpp"
#include "lattice_ops.hpp"
#include "test.hpp"

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
    // Boilerplate
    std::random_device r;
    rng.seed(r());
    lattice = new link[lsites];
    simulation1(argc,argv);



    //free memory
    // delete[] lattice;
    return 0;
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
