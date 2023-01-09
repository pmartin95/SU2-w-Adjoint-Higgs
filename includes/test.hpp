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
#pragma once

void simulation1(int argc, char** argv);
void simulation2(int argc, char** argv);
void confIdentical();
void simulation3(int argc, char** argv);