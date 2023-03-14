#include "utility.hpp"

bool verify::equivalent(const matrix &A, const matrix &B)
{
    double eps = std::numeric_limits<double>::epsilon();
    return verify::equivalent(A, B, eps);
}
bool verify::equivalent(const std::complex<double> &A, const std::complex<double> &B)
{
    return verify::equivalent(A, B, 1.0);
}
bool verify::equivalent(const double &A, const double &B)
{
    double eps = std::numeric_limits<double>::epsilon();
    return verify::equivalent(A, B, eps);
}
bool verify::equivalent(const matrix &A, const matrix &B, double epsilon)
{
    bool status = true;
    status = status && verify::equivalent(A(0, 0), B(0, 0), epsilon);
    status = status && verify::equivalent(A(0, 1), B(0, 1), epsilon);
    status = status && verify::equivalent(A(1, 0), B(1, 0), epsilon);
    status = status && verify::equivalent(A(1, 1), B(1, 1), epsilon);
    return status;
}
bool verify::equivalent(const std::complex<double> &A, const std::complex<double> &B, double epsilon)
{
    return verify::equivalent(A.real(), B.real(), epsilon) && verify::equivalent(A.imag(), B.imag(), epsilon);
}
bool verify::equivalent(const double &A, const double &B, double mult)
{
    double eps = std::numeric_limits<double>::epsilon();
    return std::abs(A - B) < eps * mult * std::max(1.0, std::max(std::abs(A), std::abs(B)));
}

void setupPauliMatrices()
{
    pauliMatrix[0] = matrix::Identity();
    pauliMatrix[1] << 0, 1, 1, 0;
    pauliMatrix[2] << 0, -I, I, 0;
    pauliMatrix[3] << 1, 0, 0, -1;
}