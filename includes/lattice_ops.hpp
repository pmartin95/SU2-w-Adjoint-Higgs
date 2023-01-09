#include "global_decl.hpp"
#pragma once
// coordinates to site index
int coordinatesToSiteIndex(int t, int x, int y, int z);
int coordinatesToSiteIndex(int (&r)[4]);
// site index to coordinates
void siteIndexToCoordinates(int site_index, int &t, int &x, int &y, int &z);
void siteIndexToCoordinates(int site_index, int (&r)[4]);
int siteJump(int site_index, int jump_dir, int jump_len, int &num_overlaps);
int siteJump(int site_index, int (&jump_len)[4], int (&num_overlaps)[4]);

void hotLattice();
void coldLattice();

// boundary condition related functions and twists
void callLatticeSite(matrix &m, int ref_site, int (&jump)[4], int dir);
matrix callLatticeSite(int ref_site, int (&jump)[4], int dir);
void ptwist(matrix &m, int mat_num, int dir, int num_twist);
void ctwist(matrix &m, int mat_num, int dir, int num_twist);
void twist(matrix &m, int mat_num, int dir, int num_twist);
void xtwist(matrix &m, int num_twist);
void ytwist(matrix &m, int num_twist);
void ztwist(matrix &m, int num_twist);