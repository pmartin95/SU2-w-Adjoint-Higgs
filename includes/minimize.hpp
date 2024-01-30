#include "simulation.hpp"
#include "action.hpp"
#include "metropolis-hastings.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#pragma once

namespace minimize{
    void simulatedAnnealingSweep(std::unique_ptr<Simulation> &sim,double T);
    void evolveLink(std::unique_ptr<Simulation> &sim,int site_index, int dir,double T);
    void evolveHiggs(std::unique_ptr<Simulation> &sim,int site_index,double T);
}