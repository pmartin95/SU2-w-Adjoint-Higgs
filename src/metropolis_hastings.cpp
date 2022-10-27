#include "global_decl.hpp"
#include "metropolis-hastings.hpp"
#include "action.hpp"
#include <iostream>

void metropolisHastingsSweep()
{
    matrix tmp;
    double actionDiff;
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < 4; dir++)
            evolveLink(site_index, dir);
        evolveHiggs(site_index);
    }
}

void evolveLink(int site_index, int dir)
{
    matrix tmp;
    double actionDiff;
    actionDiff = GeorgiGlashowPartialActionLink(site_index, dir);
    tmp = lattice[site_index].field[dir];
    generateRandomSU2Rot(lattice[site_index].field[dir]);
    actionDiff -= GeorgiGlashowPartialActionLink(site_index, dir); // actionDiff = -Delta S

    actionDiff = std::exp(actionDiff);
    if (actionDiff < gen(rng))
    {
        lattice[site_index].field[dir] = tmp;
        ++Nreject;
        ++NrejectLink;
    }
    else
    {
        ++Naccept;
        ++NacceptLink;
    }
}
void evolveHiggs(int site_index)
{
    matrix tmp;
    double actionDiff;
    actionDiff = GeorgiGlashowPartialActionHiggs(site_index);
    tmp = lattice[site_index].field[4];
    stepTracelessHermitian(lattice[site_index].field[4]);
    actionDiff -= GeorgiGlashowPartialActionHiggs(site_index); // actionDiff = -Delta S

    actionDiff = std::exp(actionDiff);
    if (actionDiff < gen(rng))
    {
        lattice[site_index].field[4] = tmp;
        ++Nreject;
        ++NrejectHiggs;
    }
    else
    {
        ++Naccept;
        ++NacceptHiggs;
    }
}