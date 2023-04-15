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
#include "utility.hpp"
#include "lattice_ops.hpp"
#pragma once

void simulation1(int argc, char **argv);
void simulation2(int argc, char **argv);
void confIdentical();
void simulation3(int argc, char **argv);
void simulation1mod(int argc, char **argv);
void simulation2mod(int argc, char **argv);
void simulation3mod(int argc, char **argv);
void simulation4(int argc, char **argv);
void argumentInput(int argc, char **argv);
int configureStep();
void boundaryConditionTest(int argc, char **argv);
void indexTest();
void simulation5(int argc, char **argv);
void midSimObservables();
void simulationHMC1(int argc, char **argv);
void testExpCK(int num_tests);