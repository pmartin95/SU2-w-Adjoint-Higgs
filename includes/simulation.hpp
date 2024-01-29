#pragma once

#include <Eigen/Dense>
#include <random>
#include <complex>
#include <string>

typedef Eigen::Matrix<std::complex<double>, 2, 2> matrix;
typedef struct site
{
    matrix field[5];
} site;

typedef void (*boundary_condition)(matrix &, int, int, int);

class Simulation
{
public:
    Simulation();  // Constructor
    ~Simulation(); // Destructor

    // Methods for operations that involve these members
    // Getters
    const std::vector<site> &getLattice() const;
    boundary_condition getBoundaryCondition() const;
    double getBeta() const;
    double getLambda() const;
    double getM2() const;
    double getKappa() const;
    int getMaxIter() const;
    int getIterCount() const;
    int getNaccept() const;
    int getNreject() const;
    int getNacceptLink() const;
    int getNrejectLink() const;
    int getNacceptHiggs() const;
    int getNrejectHiggs() const;
    int getL() const;
    int getLt() const;
    int getLsites() const;
    const int *getLdir() const;
    bool getThermalize() const;
    bool getConfigSteps() const;
    double getRotSize() const;
    double getStepSizeHiggs() const;
    const std::mt19937 &getRng() const;
    const std::uniform_real_distribution<double> &getGen() const;
    const std::normal_distribution<double> &getGenNormal() const;
    const std::string &getDatFolder() const;
    const std::string &getConfFolder() const;
    const std::string &getBcName() const;
    const std::string &getIdentifier() const;

    // Setters
    void setLattice(const std::vector<site> &lattice);
    void setBoundaryCondition(boundary_condition bc);
    void setBeta(double beta);
    void setLambda(double lambda);
    void setM2(double m2);
    void setKappa(double kappa);
    void setMaxIter(int MAX_ITER);
    void setIterCount(int iter_count);
    void setNaccept(int Naccept);
    void setNreject(int Nreject);
    void setNacceptLink(int NacceptLink);
    void setNrejectLink(int NrejectLink);
    void setNacceptHiggs(int NacceptHiggs);
    void setNrejectHiggs(int NrejectHiggs);
    void setL(int l);
    void setLt(int lt);
    void setLsites(int lsites);
    void setLdir(const int ldir[4]);
    void setThermalize(bool thermalize);
    void setConfigSteps(bool configSteps);
    void setRotSize(double rot_size);
    void setStepSizeHiggs(double step_size_higgs);
    void setRng(const std::mt19937 &rng);
    void setGen(const std::uniform_real_distribution<double> &gen);
    void setGenNormal(const std::normal_distribution<double> &gen_normal);
    void setDatFolder(const std::string &datFolder);
    void setConfFolder(const std::string &confFolder);
    void setBcName(const std::string &bcName);
    void setIdentifier(const std::string &identifier);

private:
    std::vector<site>;
    boundary_condition bc;

    double beta, lambda, m2, kappa;
    int MAX_ITER, iter_count, Naccept, Nreject, NacceptLink, NrejectLink, NacceptHiggs, NrejectHiggs;
    int l, lt, lsites, ldir[4];
    bool thermalize, configSteps;
    bool configSteps;
    const static std::complex<double> I = std::complex(0.0, 1.0);
    double rot_size, step_size_higgs;
    std::mt19937 rng;
    std::uniform_real_distribution<double> gen;
    std::normal_distribution<double> gen_normal;
    std::string datFolder, confFolder, bcName, identifier;
};

// You may also want to define the constructor and destructor where you initialize and clean up resources
Simulation::Simulation()
{
    // Initialize your variables, allocate memory, etc.
}

Simulation::~Simulation()
{
    // Clean up, deallocate memory, etc.
}
