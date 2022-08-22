// coordinates to site index
int coordinatesToSiteIndex(int t, int x, int y, int z);
// site index to coordinates
void siteIndexToCoordinates(int site_index, int &t, int &x, int &y, int &z);
// plaquettes

double action();
double actionPartial(int site_index, int mu);
void generateRandomSU2(matrix &m);
void generateRandomSU2Rot(matrix &m);

void hotLattice();
void coldLattice();
double average(const std::vector<double> &input);
double sumVec(const std::vector<double> &input);
int computeJackknifeStatistics(const std::vector<double> &inputData, int setLength, double &Jackknife_ave, double &Jackknife_error);
