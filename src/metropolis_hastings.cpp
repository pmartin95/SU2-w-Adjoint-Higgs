#include "global_decl.hpp"
#include "metropolis-hastings.hpp"
#include "generic_func.hpp"

void metropolisHastingsSweep()
{
    matrix tmp;
    double actionDiff;
    for (int site_index = 0; site_index < lsites; site_index++)
    {

        for (int dir = 0; dir < 4; dir++)
        {
            actionDiff = actionPartial(site_index, dir);
            tmp = lattice[site_index].field[dir];
            generateRandomSU2Rot(lattice[site_index].field[dir]);
            actionDiff -= actionPartial(site_index, dir); // actionDiff = -Delta S

            actionDiff = std::exp(actionDiff);
            if (actionDiff < gen(rng))
            {
                lattice[site_index].field[dir] = tmp;
                ++Nreject;
            }
            else
            {
                ++Naccept;
            }
        }
    }
}
