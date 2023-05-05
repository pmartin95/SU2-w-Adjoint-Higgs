#include "global_decl.hpp"
#include <stdbool.h>
#include <limits>
#pragma once
namespace verify
{
    bool equivalent(const matrix &A, const matrix &B, double epsilon = std::numeric_limits<double>::epsilon()); // matrix
    bool equivalent(const std::complex<double> &A, const std::complex<double> &B, double tolerance = std::numeric_limits<double>::epsilon());
    bool equivalent(const double &A, const double &B, double tolerance = std::numeric_limits<double>::epsilon());
}
bool isAntiHermitian(const matrix &m);
bool isHermitian(const matrix &m, double tolerance = std::numeric_limits<double>::epsilon());
bool isTraceless(const matrix &m);
bool isSU2(const matrix &m, double tolerance = 1.0e-6);
bool isLatticeSU2(const site *l, double tolerance = 1e-6);
void projectSU2(matrix &m, double tolerance = 1.0e-6);
void setupPauliMatrices();