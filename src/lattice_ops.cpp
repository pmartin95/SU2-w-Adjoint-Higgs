#include "lattice_ops.hpp"
#include "global_decl.hpp"
int coordinatesToSiteIndex(int t, int x, int y, int z)
{
    int site_index = t;
    site_index = site_index * l + x;
    site_index = site_index * l + y;
    site_index = site_index * l + z;

    return site_index;
}

void siteIndexToCoordinates(int site_index, int &t, int &x, int &y, int &z)
{
    int temp = site_index;
    z = temp % l;
    temp = temp / l;
    y = temp % l;
    temp = temp / l;
    x = temp % l;
    t = temp / l;
}

void hotLattice()
{
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            generateRandomSU2(lattice[site_index].field[dir]);
        }
    }
}

void coldLattice()
{
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            lattice[site_index].field[dir] = matrix::Identity();
        }
    }
}