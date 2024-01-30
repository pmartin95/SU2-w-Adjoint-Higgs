#include "simulation.hpp"
#include <memory>
#pragma once
double WilsonAction(std::unique_ptr<Simulation> &sim);
double WilsonActionPartial(std::unique_ptr<Simulation> &sim, int site_index, int mu);
double GeorgiGlashowHiggsPotential(std::unique_ptr<Simulation> &sim);
double GeorgiGlashowHiggsPotentialSite(std::unique_ptr<Simulation> &sim, int site_index);
double GeorgiGlashowHiggsKinetic(std::unique_ptr<Simulation> &sim);
double GeorgiGlashowHiggsKineticSite(std::unique_ptr<Simulation> &sim, int site_index);
double GeorgiGlashowHiggsKineticSiteOneLink(std::unique_ptr<Simulation> &sim, int site_index, int dir);
double GeorgiGlashowAction(std::unique_ptr<Simulation> &sim);
double GeorgiGlashowPartialActionLink(std::unique_ptr<Simulation> &sim, int site_index, int dir);
double GeorgiGlashowPartialActionHiggs(std::unique_ptr<Simulation> &sim, int site_index);
double HiggsMixedTerm(std::unique_ptr<Simulation>& sim,int site_index, int dir);