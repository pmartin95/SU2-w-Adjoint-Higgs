#include "simulation.hpp"
#pragma once
// coordinates to site index
int coordinatesToSiteIndex(const simulation *sim, int t, int x, int y, int z);
int coordinatesToSiteIndex(const simulation *sim, int (&r)[4]);
// site index to coordinates
void siteIndexToCoordinates(const simulation *sim, int site_index, int &t, int &x, int &y, int &z);
void siteIndexToCoordinates(const simulation *sim, int site_index, int (&r)[4]);
int siteJump(const simulation *sim, int site_index, int jump_dir, int jump_len, int &num_overlaps);
int siteJump(const simulation *sim, int site_index, int (&jump_len)[4], int (&num_overlaps)[4]);

void hotLattice(const simulation *sim);
void coldLattice(const simulation *sim);

// boundary condition related functions and twists
void callLatticeSite(const simulation *sim, matrix &m, int ref_site, int (&jump)[4], int dir);
matrix callLatticeSite(const simulation *sim, int ref_site, int (&jump)[4], int dir);
void ptwist(matrix &m, int mat_num, int dir, int num_twist);
void ctwist(matrix &m, int mat_num, int dir, int num_twist);
void twist(matrix &m, int mat_num, int dir, int num_twist);
void xtwist(matrix &m, int num_twist);
void ytwist(matrix &m, int num_twist);
void ztwist(matrix &m, int num_twist);