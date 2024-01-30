#include "rand.hpp"
#include "action.hpp"
#include "simulation.hpp"
#include <memory>
#pragma once
void metropolisHastingsSweep(std::unique_ptr<Simulation> &sim);
void evolveLink(std::unique_ptr<Simulation> &sim,int site_index, int dir);
void evolveHiggs(std::unique_ptr<Simulation> &sim,int site_index);