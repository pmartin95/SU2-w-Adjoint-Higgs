#include "hmc.hpp"

// update lattice momenta of lattice 1
void updateLatticeMomenta(site *&current_momenta, site *&current_lattice, site *&next_momenta, double step_size)
{
    for (int i = 0; i < lsites; i++)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            matrix force;
            linkFieldForce(current_lattice, i, mu, force);
            if (!isHermitian(force))
            {
                std::cout << "force is not Hermitian" << std::endl;
                exit(1);
            }
            next_momenta[i].field[mu] = current_momenta[i].field[mu] - step_size * force;
        }
    }
    std::swap(current_momenta, next_momenta);
}
// update lattice fields of lattice 1
void updateLatticeFields(site *&current_lattice, site *&current_momenta, site *&next_lattice, double step_size)
{
    for (int i = 0; i < lsites; i++)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            if (!isHermitian(current_momenta[i].field[mu]))
            {
                std::cout << "Momentum is not Hermitian" << std::endl;
                std::cout << current_momenta[i].field[mu] << std::endl;
                exit(1);
            }
            next_lattice[i].field[mu] = expCK(-step_size * I * current_momenta[i].field[mu]) * current_lattice[i].field[mu]; // Seems like  Re and + is best so far
        }
    }

    std::swap(current_lattice, next_lattice);
}

void HMC_warmup(int t)
{

    double t_step = 1.0 / t;
    site *current_lattice = lattice1_global, *next_lattice = lattice2_global;
    site *current_momenta = plattice1_global, *next_momenta = plattice2_global;

    randomMomentumLattice(plattice);
    copyLattice(lattice, current_lattice);
    copyLattice(plattice, current_momenta);

    // random momenta
    // Update lattice momenta
    updateLatticeMomenta(current_momenta, current_lattice, next_momenta, t_step / 2.0);
    // Update lattice fields
    for (int i = 0; i < t - 1; i++)
    {
        updateLatticeFields(current_lattice, current_momenta, next_lattice, t_step);
        updateLatticeMomenta(current_momenta, current_lattice, next_momenta, t_step);
    }
    updateLatticeFields(current_lattice, current_momenta, next_lattice, t_step);
    updateLatticeMomenta(current_momenta, current_lattice, next_momenta, t_step / 2.0);
    copyLattice(current_lattice, lattice);
    copyLattice(current_momenta, plattice);
}

// perform an HMC with t time intervals
double HMC(int t)
{
    double t_step = 1.0 / static_cast<double>(t);
    // Copy lattice  to lattice 1 and 2
    copyLattice(lattice, lattice1_global);
    copyLattice(lattice, lattice2_global);

    site *current_lattice = lattice1_global, *next_lattice = lattice2_global;
    site *current_momenta = plattice1_global, *next_momenta = plattice2_global;

    // random momenta
    randomMomentumLattice(current_momenta);
    copyLattice(current_momenta, plattice);

    // Update lattice momenta

    updateLatticeMomenta(current_momenta, current_lattice, next_momenta, t_step / 2.0);

    // Update lattice fields
    for (int i = 0; i < t - 1; i++)
    {
        updateLatticeFields(current_lattice, current_momenta, next_lattice, t_step);
        updateLatticeMomenta(current_momenta, current_lattice, next_momenta, t_step);
    }
    updateLatticeFields(current_lattice, current_momenta, next_lattice, t_step);
    updateLatticeMomenta(current_momenta, current_lattice, next_momenta, t_step / 2.0);
    // Calculate Hamiltonian
    double H1 = hamiltonian(lattice, plattice);
    double H2 = hamiltonian(current_lattice, current_momenta);
    std::cout << "Hdiff: " << H2 - H1 << std::endl;

    // Accept or reject
    if (H2 <= H1)
    {
        copyLattice(current_lattice, lattice);
        copyLattice(current_momenta, plattice);
        Naccept++;
    }
    else
    {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double r = distribution(rng);
        if (r < std::exp(H1 - H2))
        {
            copyLattice(current_lattice, lattice);
            copyLattice(current_momenta, plattice);
            Naccept++;
            std::cout << "Accepted exp" << std::endl;
        }
        else
        {
            Nreject++;
        }
    }
    return std::exp(H2 - H1);
}

// Copy lattice 1 to lattice 2
void copyLattice(site *lattice1, site *lattice2)
{
    for (int i = 0; i < lsites; i++)
    {
        for (int mu = 0; mu < 5; mu++)
        {
            lattice2[i].field[mu] = lattice1[i].field[mu];
        }
    }
}

// Link field force term i.e. \partial S / \partial A_{mu}
void linkFieldForce(site *lattice1, int site_index, int mu, matrix &force)
{
    int jumpNone[4] = {0, 0, 0, 0};
    matrix staple = gaugeFieldStaple(lattice1, site_index, mu);
    matrix link = callLatticeSite(lattice1, site_index, jumpNone, mu);
    matrix temp = link * staple;
    force = beta / 4.0 * antiHermitianTraceless(temp) * I;
}
const matrix linkFieldForce(site *lattice1, int site_index, int mu)
{
    matrix force;
    linkFieldForce(lattice1, site_index, mu, force);
    return force;
}

// Find the staple of the gauge field at site_index in direction mu on lattice1
void gaugeFieldStaple(site *lattice1, int site_index, int mu, matrix &staple)
{
    staple.setZero();
    int jumpNone[4] = {0, 0, 0, 0};
    for (int nu = 0; nu < 4; nu++)
    {
        if (nu != mu)
        {
            int jump1[4] = {0, 0, 0, 0};
            int jump2[4] = {0, 0, 0, 0};
            jump1[mu] = 1;
            jump2[nu] = 1;
            staple += callLatticeSite(lattice1, site_index, jump1, nu) * callLatticeSite(lattice1, site_index, jump2, mu).adjoint() * callLatticeSite(lattice1, site_index, jumpNone, nu).adjoint();
            jump1[nu] = -1;
            jump2[nu] = -1;
            staple += callLatticeSite(lattice1, site_index, jump1, nu).adjoint() * callLatticeSite(lattice1, site_index, jump2, mu).adjoint() * callLatticeSite(lattice1, site_index, jump2, nu);
        }
    }
}
const matrix gaugeFieldStaple(site *lattice1, int site_index, int mu)
{
    matrix staple;
    gaugeFieldStaple(lattice1, site_index, mu, staple);
    return staple;
}

// Cayley-Klein version of exponential function for traceless anti-Hermitian matrices
void expCK(const matrix &m, matrix &expm)
{
    std::complex<double> M = std::sqrt(m.determinant());
    expm = std::cos(M) * matrix::Identity() + std::sin(M) / M * m;
}
const matrix expCK(const matrix &m)
{
    matrix expm;
    expCK(m, expm);
    return expm;
}

double totalMomentum(site *plattice1)
{
    double total = 0.0;
    for (int i = 0; i < lsites; i++)
    {
        for (int mu = 0; mu < 5; mu++)
        {
            total += (plattice1[i].field[mu] * plattice1[i].field[mu]).trace().real();
        }
    }
    return total;
}

// Hamiltonian of the system
double hamiltonian(site *lattice1, site *plattice1)
{
    double H = 0.0;
    for (int i = 0; i < lsites; i++)
    {
        for (int mu = 0; mu < 5; mu++)
        {
            H += 0.5 * (plattice1[i].field[mu] * plattice1[i].field[mu]).trace().real();
        }
    }
    H += WilsonAction(lattice1);
    return H;
}

// generate random normally distributed Hermitian matrices, traceless
void randomHermitianMatrix(matrix &m)
{
    double a, b, c;
    std::normal_distribution<double> distribution(0.0, 1.0);
    matrix herm;
    a = distribution(rng);
    b = distribution(rng);
    c = distribution(rng);
    herm(0, 0) = std::complex<double>(c);
    herm(1, 1) = std::complex<double>(-c);
    herm(0, 1) = std::complex<double>(a, -b);
    herm(1, 0) = std::complex<double>(a, b);
    m = herm / 2.0; // Dividing by 2 set's it to be in the adjoint basis instead of the Pauli basis
}
void randomMomentumLattice(site *lattice1)
{
    for (int i = 0; i < lsites; i++)
    {
        for (int mu = 0; mu < 5; mu++)
        {
            randomHermitianMatrix(lattice1[i].field[mu]);
        }
    }
}

matrix antiHermitianTraceless(matrix &m)
{
    matrix antiherm;
    antiherm = m - m.adjoint();

    return antiherm * 0.5 - matrix::Identity() * 0.25 * (antiherm).trace();
}