#include "global_decl.hpp"
#pragma once

double WilsonAction();
double WilsonAction(site *lattice1);
double WilsonActionPartial(int site_index, int mu);
double GeorgiGlashowHiggsPotential();
double GeorgiGlashowHiggsPotentialSite(int site_index);
double GeorgiGlashowHiggsKinetic();
double GeorgiGlashowHiggsKineticSite(int site_index);
double GeorgiGlashowHiggsKineticSiteOneLink(int site_index, int dir);
double GeorgiGlashowAction();
double GeorgiGlashowPartialActionLink(int site_index, int dir);
double GeorgiGlashowPartialActionHiggs(int site_index);
double HiggsMixedTerm(int site_index, int dir);