#include "simulation.hpp"
#include <memory>
#pragma once
// coordinates to site index
int coordinatesToSiteIndex(std::unique_ptr<Simulation> &sim, int t, int x, int y, int z);
int coordinatesToSiteIndex(std::unique_ptr<Simulation> &sim, int (&r)[4]);
// site index to coordinates
void siteIndexToCoordinates(std::unique_ptr<Simulation> &sim, int site_index, int &t, int &x, int &y, int &z);
void siteIndexToCoordinates(std::unique_ptr<Simulation> &sim, int site_index, int (&r)[4]);
int siteJump(std::unique_ptr<Simulation> &sim, int site_index, int jump_dir, int jump_len, int &num_overlaps);
int siteJump(std::unique_ptr<Simulation> &sim, int site_index, int (&jump_len)[4], int (&num_overlaps)[4]);

void hotLattice(std::unique_ptr<Simulation> &sim);
void coldLattice(std::unique_ptr<Simulation> &sim);

// boundary condition related functions and twists
void callLatticeSite(std::unique_ptr<Simulation> &sim, matrix &m, int ref_site, int (&jump)[4], int dir);
matrix callLatticeSite(std::unique_ptr<Simulation> &sim, int ref_site, int (&jump)[4], int dir);
void ptwist(matrix &m, int mat_num, int dir, int num_twist);
void ctwist(matrix &m, int mat_num, int dir, int num_twist);
void twist(matrix &m, int mat_num, int dir, int num_twist);
void xtwist(matrix &m, int num_twist);
void ytwist(matrix &m, int num_twist);
void ztwist(matrix &m, int num_twist);