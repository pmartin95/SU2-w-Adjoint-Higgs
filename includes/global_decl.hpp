#pragma once
#include <Eigen/Dense>
#include <random>
#include <complex>
#include <string>
typedef Eigen::Matrix<std::complex<double>, 2, 2> matrix;
typedef struct link
{
    matrix field[5];
} link;
extern link *lattice;

extern double beta;
extern double lambda;
extern double m2;
extern double kappa;
extern int MAX_ITER;
extern int iter_count;
extern int Naccept;
extern int Nreject;
extern int NacceptLink;
extern int NrejectLink;
extern int NacceptHiggs;
extern int NrejectHiggs;
extern int l;
extern int lt;
extern int lsites;
extern int ldir[4];
extern std::complex<double> I;
extern double rot_size;
extern double step_size_higgs;
extern std::mt19937 rng;
extern std::uniform_real_distribution<double> gen;
extern std::normal_distribution<double> gen_normal;
extern std::string datFolder;
extern std::string confFolder;
