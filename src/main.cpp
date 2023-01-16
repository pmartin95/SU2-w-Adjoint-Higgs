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

#include <boost/filesystem.hpp>
#include <boost/version.hpp>

#include "global_decl.hpp"
#include "observables.hpp"
#include "metropolis-hastings.hpp"
#include "statistics.hpp"
#include "configuration_io.hpp"
#include "lattice_ops.hpp"
#include "test.hpp"
#include "rand.hpp"
#include "action.hpp"
site *lattice;
double beta = 2.3;
double lambda = 0.1;
double m2 = 2.0;
double kappa = 0.0; // This term controls the strength of the Higgs kinetic term in the action.
int MAX_ITER = 10;
int iter_count = 0;
int Naccept = 0, Nreject = 0;
int NacceptLink = 0, NrejectLink = 0;
int NacceptHiggs = 0, NrejectHiggs = 0;
int l = 4, lt = 4, lsites = l * l * l * lt;
int ldir[4] = {lt, l, l, l};
bool thermalize = true;
bool configSteps = true;
std::complex<double> I(0.0, 1.0);
double rot_size = 0.4;
double step_size_higgs = 0.07;
std::mt19937 rng;
std::uniform_real_distribution<double> gen(0.0, 1.0);
std::normal_distribution<double> gen_normal(0.0, 1.0);
boundary_condition bc = &ptwist;
std::string datFolder = "./dat/testdata6/";
std::string confFolder = "./configurations/";
std::string bcName = "p";
std::string identifier;
int main(int argc, char **argv)
{
    if (!boost::filesystem::exists(boost::filesystem::path(datFolder)))
        boost::filesystem::create_directory(boost::filesystem::path(datFolder));
    // Boilerplate
    std::random_device r;
    rng.seed(r());
    lattice = new site[lsites];
    simulation1(argc, argv);
    delete[] lattice;
    return 0;
}
