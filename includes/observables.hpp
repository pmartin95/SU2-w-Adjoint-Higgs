#include <string>
#include "global_decl.hpp"
#pragma once
double plaquette(int site_index, int (&shift)[4], int mu, int nu);
double plaquetteAverage();
double polyakovLine(int site_index);
double polyakovLine(int site_index, int dir);
void polyakovLines(std::string filename, int dir1, int dir2);
void polyakovLinesAbs(std::string filename);
double rectangle(int site_index, int mu, int nu, int mu_len, int nu_len);
double rectangleAverage(int mu_len, int nu_len);
double higgsSquareAverage();
matrix higgsAverage();
double correlator(int site_index, int time_forward);
double averageCorrelatorVolume(int time_forward);

void copyCoordinates(int (&src)[4], int (&dest)[4]);
void addCoordinatesInPlace(int (&src)[4], int (&dest)[4]);