#include "configuration_io.hpp"
#include "global_decl.hpp"
#include <fstream>
const int MATRICES_PER_SITE = 5;
const int DOUBLES_PER_MATRIX = 8;

int pushConfig(std::string filename)
{
    std::vector<double> v;
    v.reserve(lsites * MATRICES_PER_SITE * DOUBLES_PER_MATRIX);
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < MATRICES_PER_SITE; dir++)
        {
            v.push_back(lattice[site_index].field[dir](0, 0).real());
            v.push_back(lattice[site_index].field[dir](0, 0).imag());
            v.push_back(lattice[site_index].field[dir](0, 1).real());
            v.push_back(lattice[site_index].field[dir](0, 1).imag());
            v.push_back(lattice[site_index].field[dir](1, 0).real());
            v.push_back(lattice[site_index].field[dir](1, 0).imag());
            v.push_back(lattice[site_index].field[dir](1, 1).real());
            v.push_back(lattice[site_index].field[dir](1, 1).imag());
        }
    }

    std::ofstream conf(filename, std::ios::binary);
    if (!conf)
    {
        return -1;
    }
    // conf.write(reinterpret_cast<const char *>(&v[0]), v.size() * sizeof(v[0])); old way of doing it
    conf.write(reinterpret_cast<const char *>(v.data()), v.size() * sizeof(v[0]));
    return 0;
}
int pullConfig(std::string filename)
{
    int i, j, k;
    int doublesPerSite = MATRICES_PER_SITE * DOUBLES_PER_MATRIX;

    std::vector<double> v(lsites * doublesPerSite);
    std::ifstream conf;
    conf.open(filename, std::ios::binary);
    if (!conf)
    {
        return -1;
    }
    // conf.read(reinterpret_cast<char *>(&v[0]), lsites * doublesPerSite * sizeof(double)); // old way
    conf.read(reinterpret_cast<char *>(v.data()), lsites * doublesPerSite * sizeof(double));

    for (int site_index = 0; site_index < lsites; site_index++)
    {
        k = site_index * doublesPerSite;
        j = 0;
        for (int i = 0; i < 5; i++)
        {

            lattice[site_index].field[i](0, 0) = std::complex<double>(v[k + j], v[k + j + 1]);
            j += 2;
            lattice[site_index].field[i](0, 1) = std::complex<double>(v[k + j], v[k + j + 1]);
            j += 2;
            lattice[site_index].field[i](1, 0) = std::complex<double>(v[k + j], v[k + j + 1]);
            j += 2;
            lattice[site_index].field[i](1, 1) = std::complex<double>(v[k + j], v[k + j + 1]);
            j += 2;
        }
    }
    return 0;
}

void createIdentifier(std::string &unique_key)
{
    unique_key = "m2_" + std::to_string(m2) + "beta" + std::to_string(beta) +
                 "lambda" + std::to_string(lambda) + "kappa" + std::to_string(kappa) +
                 "l" + std::to_string(l) + "t" + std::to_string(lt) + "bc" + bcName;
}