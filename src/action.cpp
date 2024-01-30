#include "action.hpp"
#include "Simulation.hpp"
#include "observables.hpp"
#include "lattice_ops.hpp"
#include <memory>

double WilsonAction(std::unique_ptr<Simulation> &sim)
{
    double accumulator = 0.0;
    int jumpNone[4] = {0};
    for (int site_index = 0; site_index < sim->lsites; site_index++)
    {
        for (int nu = 0; nu < 4; nu++)
        {
            for (int mu = 0; mu < nu; mu++)
            {
                accumulator += (1.0 - plaquette(sim, site_index, jumpNone, mu, nu));
            }
        }
    }
    return sim->beta * accumulator;
}
double WilsonActionPartial(std::unique_ptr<Simulation> &sim, int site_index, int mu)
{
    double accumulator = 0.0;
    int jumpNone[4] = {0};
    for (int nu = 0; nu < 4; nu++)
    {
        if (mu == nu)
            continue;
        int x[4];
        siteIndexToCoordinates(sim, site_index, x[0], x[1], x[2], x[3]);
        x[nu] = (x[nu] - 1 + sim->ldir[nu]) % sim->ldir[nu];
        accumulator += plaquette(sim, site_index, jumpNone, mu, nu);
        accumulator += plaquette(sim, coordinatesToSiteIndex(sim, x[0], x[1], x[2], x[3]), jumpNone, mu, nu);
    }
    return -sim->beta * accumulator;
}

double GeorgiGlashowHiggsPotential(std::unique_ptr<Simulation> &sim)
{
    double accumulatorLambda = 0.0;
    double accumulatorm2 = 0.0;
    double higgsSquare;
    matrix tmp;

    for (int site_index = 0; site_index < sim->lsites; site_index++)
    {
        tmp = sim->lattice[site_index].field[4];
        higgsSquare = (tmp * tmp).trace().real();
        accumulatorLambda += higgsSquare * higgsSquare;
        accumulatorm2 += higgsSquare;
    }
    return sim->m2 * accumulatorm2 + sim->lambda * accumulatorLambda;
}

double GeorgiGlashowHiggsPotentialSite(std::unique_ptr<Simulation> &sim, int site_index)
{
    double tmp = (sim->lattice[site_index].field[4] * sim->lattice[site_index].field[4]).trace().real();
    return tmp * tmp * sim->lambda + tmp * sim->m2;
}

//! Investigate this term
double GeorgiGlashowHiggsKinetic(std::unique_ptr<Simulation> &sim)
{
    double accumulator = 0.0;
    matrix tmp;
    for (int site_index = 0; site_index < sim->lsites; site_index++)
    {
        tmp = sim->lattice[site_index].field[4];
        accumulator += 4.0 * (tmp * tmp).trace().real();
        for (int dir = 0; dir < 4; dir++)
        {
            accumulator -= HiggsMixedTerm(sim, site_index, dir);
        }
    }
    return 2.0 * sim->kappa * accumulator;
}
//! Investigate this term
double GeorgiGlashowHiggsKineticSite(std::unique_ptr<Simulation> &sim, int site_index)
{
    double accumulator;
    matrix tmp = sim->lattice[site_index].field[4];
    accumulator = 4.0 * (tmp * tmp).trace().real();

    for (int dir = 0; dir < 4; dir++)
    {
        int x[4];
        siteIndexToCoordinates(sim, site_index, x[0], x[1], x[2], x[3]);
        x[dir] = (x[dir] - 1 + sim->ldir[dir]) % sim->ldir[dir]; //! Wait a seconsd, should this be here?
        accumulator -= HiggsMixedTerm(sim, site_index, dir);
        accumulator -= HiggsMixedTerm(sim, coordinatesToSiteIndex(sim, x[0], x[1], x[2], x[3]), dir);
    }
    return 2.0 * sim->kappa * accumulator;
}
//! Investigate this term
double GeorgiGlashowHiggsKineticSiteOneLink(std::unique_ptr<Simulation> &sim, int site_index, int dir)
{
    return -2.0 * sim->kappa * HiggsMixedTerm(sim, site_index, dir);
}
double GeorgiGlashowAction(std::unique_ptr<Simulation> &sim)
{
    return WilsonAction(sim) + GeorgiGlashowHiggsPotential(sim) + GeorgiGlashowHiggsKinetic(sim);
}

double GeorgiGlashowPartialActionLink(std::unique_ptr<Simulation> &sim, int site_index, int dir)
{
    return WilsonActionPartial(sim, site_index, dir) + GeorgiGlashowHiggsKineticSiteOneLink(sim, site_index, dir);
}
double GeorgiGlashowPartialActionHiggs(std::unique_ptr<Simulation> &sim, int site_index)
{
    return GeorgiGlashowHiggsKineticSite(sim, site_index) + GeorgiGlashowHiggsPotentialSite(sim, site_index);
}
//! Investigate this term
double HiggsMixedTerm(std::unique_ptr<Simulation> &sim, int site_index, int dir) // dir = mu,  Tr( \Phi(x) U_{\mu}(x)  \Phi(x+\mu)  U_{\mu}^{\adjoint}(x) )
{
    int x[4], jump_index;

    siteIndexToCoordinates(sim, site_index, x[0], x[1], x[2], x[3]);
    x[dir] = (x[dir] + 1) % sim->ldir[dir];
    jump_index = coordinatesToSiteIndex(sim, x[0], x[1], x[2], x[3]);
    matrix higgsX, higgsXpMU, link;
    higgsX = sim->lattice[site_index].field[4];
    link = sim->lattice[site_index].field[dir];
    higgsXpMU = sim->lattice[jump_index].field[4];
    return (higgsX * link * higgsXpMU * link.adjoint()).trace().real();
}