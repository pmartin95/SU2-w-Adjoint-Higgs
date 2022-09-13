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
#include "metropolis-hastings.hpp"
#include "statistics.hpp"
#include "configuration_io.hpp"
#include "lattice_ops.hpp"
#include "test.hpp"
#include "rand.hpp"
#include "action.hpp"
link *lattice;
double beta = 2.3;
double lambda = 0.1;
double m2 = -.2;
int MAX_ITER = 1;
int iter_count = 0;
int Naccept = 0, Nreject = 0;
int NacceptLink = 0, NrejectLink = 0;
int NacceptHiggs = 0, NrejectHiggs = 0;
int l = 2, lt = 2, lsites = l * l * l * lt;
int ldir[4] = {lt, l, l, l};
std::complex<double> I(0.0, 1.0);
double rot_size = 0.4;
double step_size_higgs = 0.4;
std::mt19937 rng;
std::uniform_real_distribution<double> gen(0.0, 1.0);
std::normal_distribution<double> gen_normal(0.0, 1.0);

int main(int argc, char **argv)
{
    // Boilerplate
    std::random_device r;
    rng.seed(r());

    // // simulation1(argc,argv);
    lattice = new link[lsites];
    step_size_higgs = 0.07;
    hotLattice();
    beta = 5;
    for (int i = 0; i < 1000; i++)
    {
        metropolisHastingsSweep();
        std::cout << higgsSquareAverage() << std::endl;
    }
    delete[] lattice;
    // // free memory

    return 0;
}
