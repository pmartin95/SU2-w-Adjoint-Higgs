#include "minimize.hpp"
namespace minimize
{
    void simulatedAnnealingSweep(std::unique_ptr<Simulation> &sim, double T)
    {
        double actionDiff;
        for (int site_index = 0; site_index < sim->lsites; site_index++)
        {
            for (int dir = 0; dir < 4; dir++)
                minimize::evolveLink(sim, site_index, dir, T);
            minimize::evolveHiggs(sim, site_index, T);
        }
    }
    void evolveLink(std::unique_ptr<Simulation> &sim, int site_index, int dir, double T)
    {
        matrix tmp;
        double actionDiff;
        actionDiff = GeorgiGlashowPartialActionLink(sim, site_index, dir);
        tmp = sim->lattice[site_index].field[dir];
        generateRandomSU2Rot(sim->lattice[site_index].field[dir]);
        actionDiff -= GeorgiGlashowPartialActionLink(sim, site_index, dir); // actionDiff = -Delta S

        actionDiff = std::exp(actionDiff / T);
        if (actionDiff < sim->gen(sim->rng))
        {
            sim->lattice[site_index].field[dir] = tmp;
            sim->Nreject++;
            sim->NrejectLink++;
        }
        else
        {
            sim->Naccept++;
            sim->NacceptLink++;
        }
    }
    void evolveHiggs(std::unique_ptr<Simulation> &sim,int site_index, double T)
    {
        matrix tmp;
        double actionDiff;
        actionDiff = GeorgiGlashowPartialActionHiggs(sim, site_index);
        tmp = sim->lattice[site_index].field[4];
        stepTracelessHermitian(sim->lattice[site_index].field[4]);
        actionDiff -= GeorgiGlashowPartialActionHiggs(sim, site_index); // actionDiff = -Delta S

        actionDiff = std::exp(actionDiff / T);
        if (actionDiff < sim->gen(sim->rng))
        {
            sim->lattice[site_index].field[4] = tmp;
            sim->Nreject++;
            sim->NrejectHiggs++;
        }
        else
        {
            sim->Naccept++;
            sim->NacceptHiggs++;
        }
    }
}
