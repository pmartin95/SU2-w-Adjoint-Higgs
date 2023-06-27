#include "global_decl.hpp"
#include "action.hpp"
#include "metropolis-hastings.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#pragma once

namespace minimize{
    void simulatedAnnealingSweep(double T);
    void evolveLink(int site_index, int dir,double T);
    void evolveHiggs(int site_index,double T);
}