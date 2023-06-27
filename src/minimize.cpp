#include "minimize.hpp"
namespace minimize{
void simulatedAnnealingSweep(double T)
{
    double actionDiff;
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < 4; dir++)
            minimize::evolveLink(site_index, dir,T);
        minimize::evolveHiggs(site_index,T);


    }
}
void evolveLink(int site_index, int dir,double T)
{
    matrix tmp;
    double actionDiff;
    actionDiff = GeorgiGlashowPartialActionLink(site_index, dir);
    tmp = lattice[site_index].field[dir];
    generateRandomSU2Rot(lattice[site_index].field[dir]);
    actionDiff -= GeorgiGlashowPartialActionLink(site_index, dir); // actionDiff = -Delta S

    actionDiff = std::exp(actionDiff/T);
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
void evolveHiggs(int site_index,double T)
{
    matrix tmp;
    double actionDiff;
    actionDiff = GeorgiGlashowPartialActionHiggs(site_index);
    tmp = lattice[site_index].field[4];
    stepTracelessHermitian(lattice[site_index].field[4]);
    actionDiff -= GeorgiGlashowPartialActionHiggs(site_index); // actionDiff = -Delta S

    actionDiff = std::exp(actionDiff/T);
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
}

