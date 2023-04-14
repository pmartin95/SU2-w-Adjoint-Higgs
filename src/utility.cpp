#include "utility.hpp"
#include <iostream>

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

// Function to check if a matrix is an SU(2) matrix
bool isSU2(const matrix &m, double tolerance)
{
    // Check if the determinant is close to 1
    if (std::abs(m.determinant().real() - 1.0) > tolerance)
    {
        std::cout << m << std::endl;
        std::cout << "Determinant is not close to 1: " << m.determinant() - 1.0 << std::endl;
        return false;
    }

    // Check if the matrix is unitary: m * m.adjoint() should be the identity matrix
    matrix product = m * m.adjoint();
    if (!verify::equivalent(product, matrix::Identity(), tolerance / std::numeric_limits<double>::epsilon()))
    {
        std::cout << "Matrix is not unitary: " << std::endl;
        std::cout << product - matrix::Identity() << std::endl;
        return false;
    }

    return true;
}
bool isLatticeSU2(const site *l, double tolerance)
{
    for (int i = 0; i < lsites; i++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            if (!isSU2(l[i].field[dir], tolerance))
            {
                return false;
            }
        }
    }
    return true;
}

void projectSU2(matrix &m, double tolerance)
{
    // Project the matrix onto SU(2)
    std::complex<double> det = m.determinant();
    while (std::abs(det.real() - 1.0) > tolerance && det.real() > 0.0)
    {
        m /= std::sqrt(det);
    }
}