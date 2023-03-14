#include "global_decl.hpp"
#include <stdbool.h>
#include <limits>
#pragma once
namespace verify
{
    bool equivalent(const matrix &A, const matrix &B); // matrix
    bool equivalent(const std::complex<double> &A, const std::complex<double> &B);
    bool equivalent(const double &A, const double &B);
    bool equivalent(const matrix &A, const matrix &B, double epsilon); // matrix
    bool equivalent(const std::complex<double> &A, const std::complex<double> &B, double epsilon);
    bool equivalent(const double &A, const double &B, double mult);
}

void setupPauliMatrices();