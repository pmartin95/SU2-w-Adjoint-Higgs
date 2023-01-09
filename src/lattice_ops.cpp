#include "lattice_ops.hpp"
#include "global_decl.hpp"
#include "rand.hpp"
#include <iostream>
#include <complex>

int coordinatesToSiteIndex(int t, int x, int y, int z)
{
    int site_index = t;
    site_index = site_index * l + x;
    site_index = site_index * l + y;
    site_index = site_index * l + z;

    return site_index;
}

int coordinatesToSiteIndex(int (&r)[4])
{
    return coordinatesToSiteIndex(r[0], r[1], r[2], r[3]);
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

void siteIndexToCoordinates(int site_index, int (&r)[4])
{
    siteIndexToCoordinates(site_index, r[0], r[1], r[2], r[3]);
}

int siteJump(int site_index, int jump_dir, int jump_len, int &num_overlaps) // Returns new site index
{
    int r[4];
    num_overlaps = 0;
    siteIndexToCoordinates(site_index, r);
    r[jump_dir] += jump_len;
    while (r[jump_dir] >= ldir[jump_dir])
    {
        r[jump_dir] -= ldir[jump_dir];
        num_overlaps++;
    }
    while (r[jump_dir] < 0)
    {
        r[jump_dir] += ldir[jump_dir];
        num_overlaps--;
    }
    return coordinatesToSiteIndex(r);
}

int siteJump(int site_index, int (&jump_len)[4], int (&num_overlaps)[4]) // Returns new site index
{
    int r[4];

    siteIndexToCoordinates(site_index, r);
    for (int dir = 0; dir < 4; dir++)
    {
        num_overlaps[dir] = 0;
        r[dir] += jump_len[dir];
        while (r[dir] >= ldir[dir])
        {
            r[dir] -= ldir[dir];
            num_overlaps[dir]++;
        }
        while (r[dir] < 0)
        {
            r[dir] += ldir[dir];
            num_overlaps[dir]--;
        }
    }

    return coordinatesToSiteIndex(r);
}
void hotLattice()
{
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            generateRandomSU2(lattice[site_index].field[dir]);
        }
        // generateRandomTracelessHermitian(lattice[site_index].field[4]);
        // lattice[site_index].field[4] *= (std::abs(m2) / (2.0 * lambda));
        lattice[site_index].field[4] = matrix::Zero();
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
        lattice[site_index].field[4] << (m2 / (2.0 * lambda)), 0.0, 0.0, -(m2 / (2.0 * lambda));
    }
}

// You should rewrite variable names to make it more apparent which is being used, dir, mat_num, etc
//! call by reference function to pull field from lattice site and do boundary conditions
void callLatticeSite(matrix &m, int ref_site, int (&jump)[4], int dir)
{
    int num_twists[4];
    int new_site = siteJump(ref_site, jump, num_twists);
    m = lattice[new_site].field[dir];
    for (int i = 0; i < 4; i++)
    {
        bc(m, dir, i, num_twists[i]);
    }
}

// return copy of field from lattice site and do boundary conditions
matrix callLatticeSite(int ref_site, int (&jump)[4], int dir)
{
    matrix m;
    callLatticeSite(m, ref_site, jump, dir);
    return m;
}
// dir refers to the direction traveled
void ptwist(matrix &m, int mat_num, int dir, int num_twist) // aka no twist, just periodic b.c.
{
}
void ctwist(matrix &m, int mat_num, int dir, int num_twist)
{
    ytwist(m, num_twist);
    if (mat_num == 4 && dir != 0)
    {
        m = -m;
    }
}

void twist(matrix &m, int mat_num, int dir, int num_twist)
{
    if (dir == 1)
    {
        xtwist(m, num_twist);
    }
    else if (dir == 2)
    {
        ytwist(m, num_twist);
    }
    else if (dir == 3)
    {
        ztwist(m, num_twist);
    }
    if (mat_num == 4 && dir != 0)
    {
        m = -m;
    }
}

// \sigma_1 * m * \sigma_1
void xtwist(matrix &m, int num_twist)
{
    if (num_twist % 2 == 1)
    {
        std::complex<double> temp = m(0, 0);
        m(0, 0) = m(1, 1);
        m(1, 1) = temp;

        temp = m(0, 1);
        m(0, 1) = m(1, 0);
        m(1, 0) = temp;
    }
}
// \sigma_2 * m * \sigma_2
void ytwist(matrix &m, int num_twist)
{
    if (num_twist % 2 == 1)
    {
        std::complex<double> temp = m(0, 0);
        m(0, 0) = m(1, 1);
        m(1, 1) = temp;
        temp = -m(0, 1);
        m(0, 1) = -m(1, 0);
        m(1, 0) = temp;
    }
}
// \sigma_3 * m * \sigma_3
void ztwist(matrix &m, int num_twist)
{
    if (num_twist % 2 == 1)
    {
        m(0, 1) = -m(0, 1);
        m(1, 0) = -m(1, 0);
    }
}