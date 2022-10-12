#include "action.hpp"
#include "observables.hpp"
#include "lattice_ops.hpp"

double WilsonAction()
{
    double accumulator = 0.0;
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            for (int mu = 0; mu < nu; mu++)
            {
                accumulator += (1.0 - plaquette(site_index, mu, nu));
            }
        }
    }
    return beta * accumulator;
}
double WilsonActionPartial(int site_index, int mu)
{
    double accumulator = 0.0;
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu)
            continue;
        int x[4];
        siteIndexToCoordinates(site_index, x[0], x[1], x[2], x[3]);
        x[nu] = (x[nu] - 1 + ldir[nu]) % ldir[nu];
        accumulator += plaquette(site_index, mu, nu);
        accumulator += plaquette(coordinatesToSiteIndex(x[0], x[1], x[2], x[3]), mu, nu);
    }
    return -beta * accumulator;
}

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

double GeorgiGlashowHiggsPotentialSite(int site_index)
{
    double tmp = (lattice[site_index].field[4] * lattice[site_index].field[4]).trace().real();
    return tmp * tmp * lambda + tmp * m2;
}

//! Investigate this term
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
double GeorgiGlashowHiggsKineticSiteOneLink(int site_index, int dir)
{
    return -2.0 * kappa * HiggsMixedTerm(site_index, dir);
}
double GeorgiGlashowAction()
{
    return WilsonAction() + GeorgiGlashowHiggsPotential() + GeorgiGlashowHiggsKinetic();
}

double GeorgiGlashowPartialActionLink(int site_index, int dir)
{
    return WilsonActionPartial(site_index, dir) + GeorgiGlashowHiggsKineticSiteOneLink(site_index, dir);
}
double GeorgiGlashowPartialActionHiggs(int site_index)
{
    return GeorgiGlashowHiggsKineticSite(site_index) + GeorgiGlashowHiggsPotentialSite(site_index);
}
//! Investigate this term
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