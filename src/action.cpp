#include "action.hpp"
#include "observables.hpp"
#include "utility.hpp"
#include "lattice_ops.hpp"
#include <iostream>

// Calculate the Wilson action of the current lattice configuration
double WilsonAction()
{
    double accumulator = 0.0;
    int jumpNone[4] = {0};
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            for (int mu = 0; mu < nu; mu++)
            {
                accumulator += (1.0 - plaquette(site_index, jumpNone, mu, nu));
            }
        }
    }
    return beta * accumulator;
}

// Calculate the Wilson action of a given lattice configuration
double WilsonAction(site *lattice1)
{
    double accumulator = 0.0;
    int jumpNone[4] = {0};
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            for (int mu = 0; mu < nu; mu++)
            {
                if (std::abs(1.0 - plaquette(lattice1, site_index, jumpNone, mu, nu)) > 2.0)
                {
                    std::cout << "Plaquette is behaving strangely.\n";
                    std::cout << 2.0 - plaquette(lattice1, site_index, jumpNone, mu, nu) << "\n";
                    std::cout << isLatticeSU2(lattice1) << "\n";
                    exit(1);
                }

                accumulator += (1.0 - plaquette(lattice1, site_index, jumpNone, mu, nu));
            }
        }
    }
    return beta * accumulator;
}

// Calculate the partial Wilson action for a specific site and direction
double WilsonActionPartial(int site_index, int mu)
{
    double accumulator = 0.0;
    int jumpNone[4] = {0};
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu)
            continue;
        int x[4];
        siteIndexToCoordinates(site_index, x[0], x[1], x[2], x[3]);
        x[nu] = (x[nu] - 1 + ldir[nu]) % ldir[nu];
        accumulator += plaquette(site_index, jumpNone, mu, nu);
        accumulator += plaquette(coordinatesToSiteIndex(x[0], x[1], x[2], x[3]), jumpNone, mu, nu);
    }
    return -beta * accumulator;
}

// Calculate the Georgi-Glashow Higgs potential for the current lattice configuration
double GeorgiGlashowHiggsPotential()
{
    double accumulatorLambda = 0.0;
    double accumulatorm2 = 0.0;
    double higgsSquare;
    matrix tmp;

    for (int site_index = 0; site_index < lsites; site_index++)
    {
        tmp = lattice[site_index].field[4];
        higgsSquare = (tmp * tmp).trace().real();
        accumulatorLambda += higgsSquare * higgsSquare;
        accumulatorm2 += higgsSquare;
    }
    return m2 * accumulatorm2 + lambda * accumulatorLambda;
}

// Calculate the Georgi-Glashow Higgs potential for a specific site
double GeorgiGlashowHiggsPotentialSite(int site_index)
{
    double tmp = (lattice[site_index].field[4] * lattice[site_index].field[4]).trace().real();
    return tmp * tmp * lambda + tmp * m2;
}

//! Investigate this term
// Calculate the Georgi-Glashow Higgs kinetic term for the current lattice configuration
double GeorgiGlashowHiggsKinetic()
{
    double accumulator = 0.0;
    matrix tmp;
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        tmp = lattice[site_index].field[4];
        accumulator += 4.0 * (tmp * tmp).trace().real();
        for (int dir = 0; dir < 4; dir++)
        {
            accumulator -= HiggsMixedTerm(site_index, dir);
        }
    }
    return 2.0 * kappa * accumulator;
}

//! Investigate this term
// Calculate the Georgi-Glashow Higgs kinetic term for a specific site
double GeorgiGlashowHiggsKineticSite(int site_index)
{
    double accumulator;
    matrix tmp = lattice[site_index].field[4];
    accumulator = 4.0 * (tmp * tmp).trace().real();

    for (int dir = 0; dir < 4; dir++)
    {
        int x[4];
        siteIndexToCoordinates(site_index, x[0], x[1], x[2], x[3]);
        x[dir] = (x[dir] - 1 + ldir[dir]) % ldir[dir];
        accumulator -= HiggsMixedTerm(site_index, dir);
        accumulator -= HiggsMixedTerm(coordinatesToSiteIndex(x[0], x[1], x[2], x[3]), dir);
    }
    return 2.0 * kappa * accumulator;
}

//! Investigate this term
// Calculate the Georgi-Glashow Higgs kinetic term for a specific site and link direction
double GeorgiGlashowHiggsKineticSiteOneLink(int site_index, int dir)
{
    return -2.0 * kappa * HiggsMixedTerm(site_index, dir);
}

// Calculate the total Georgi-Glashow action for the current lattice configuration
double GeorgiGlashowAction()
{
    return WilsonAction() + GeorgiGlashowHiggsPotential() + GeorgiGlashowHiggsKinetic();
}

// Calculate the partial Georgi-Glashow action for a specific site and link direction
double GeorgiGlashowPartialActionLink(int site_index, int dir)
{
    return WilsonActionPartial(site_index, dir) + GeorgiGlashowHiggsKineticSiteOneLink(site_index, dir);
}

// Calculate the partial Georgi-Glashow action for a specific site
double GeorgiGlashowPartialActionHiggs(int site_index)
{
    return GeorgiGlashowHiggsKineticSite(site_index) + GeorgiGlashowHiggsPotentialSite(site_index);
}

//! Investigate this term
// Calculate the Higgs mixed term for a specific site and direction
double HiggsMixedTerm(int site_index, int dir) // dir = mu,  Tr( \Phi(x) U_{\mu}(x)  \Phi(x+\mu)  U_{\mu}^{\adjoint}(x) )
{
    int x[4], jump_index;

    siteIndexToCoordinates(site_index, x[0], x[1], x[2], x[3]);
    x[dir] = (x[dir] + 1) % ldir[dir];
    jump_index = coordinatesToSiteIndex(x[0], x[1], x[2], x[3]);
    matrix higgsX, higgsXpMU, link;
    higgsX = lattice[site_index].field[4];
    link = lattice[site_index].field[dir];
    higgsXpMU = lattice[jump_index].field[4];
    return (higgsX * link * higgsXpMU * link.adjoint()).trace().real();
}