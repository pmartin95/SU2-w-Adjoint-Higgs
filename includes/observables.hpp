#include <string>
#include <memory>
#include "simulation.hpp"

#pragma once
double plaquette(std::unique_ptr<Simulation> &sim, int site_index, int (&shift)[4], int mu, int nu);
double plaquetteAverage(std::unique_ptr<Simulation> &sim);
double polyakovLine(std::unique_ptr<Simulation> &sim, int site_index);
double polyakovLine(std::unique_ptr<Simulation> &sim, int site_index, int dir);
void polyakovLines(std::unique_ptr<Simulation> &sim, std::string filename, int dir1, int dir2);
void polyakovLinesAbs(std::unique_ptr<Simulation> &sim, std::string filename);
double rectangle(std::unique_ptr<Simulation> &sim, int site_index, int (&shift)[4], int mu, int nu, int mu_len, int nu_len);
double rectangleAverage(std::unique_ptr<Simulation> &sim, int mu_len, int nu_len);
double higgsSquareAverage(std::unique_ptr<Simulation> &sim);
double higgsSquareSum(std::unique_ptr<Simulation> &sim);
matrix higgsAverage(std::unique_ptr<Simulation> &sim);
double correlator(std::unique_ptr<Simulation> &sim, int site_index, int time_forward);
double averageCorrelatorVolume(std::unique_ptr<Simulation> &sim, int time_forward);

void copyCoordinates(int (&src)[4], int (&dest)[4]);
void addCoordinatesInPlace(int (&src)[4], int (&dest)[4]);