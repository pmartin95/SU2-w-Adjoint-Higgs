#include "configuration_io.hpp"
#include "global_decl.hpp"
#include <fstream>

int pushConfig(std::string filename)
{
    std::vector<double> v;
    std::ofstream conf;
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < 5; dir++)
        {
            v.push_back(lattice[site_index].field[dir](0, 0).imag());
            v.push_back(lattice[site_index].field[dir](0, 0).real());
            v.push_back(lattice[site_index].field[dir](0, 1).imag());
            v.push_back(lattice[site_index].field[dir](0, 1).real());
            v.push_back(lattice[site_index].field[dir](1, 0).imag());
            v.push_back(lattice[site_index].field[dir](1, 0).real());
            v.push_back(lattice[site_index].field[dir](1, 1).imag());
            v.push_back(lattice[site_index].field[dir](1, 1).real());
        }
    }
    conf.open(filename, std::ios::binary);
    conf.write(reinterpret_cast<const char *>(&v[0]), v.size() * sizeof(v[0]));
    conf.close();
    return 0;
}
int pullConfig(std::string filename)
{
    int i, j, k;
    int doublesPerSite = 5 /* matrices per site */ * 8 /*doubles per matrix*/;
    std::ifstream conf;
    std::vector<double> v(lsites * doublesPerSite);

    conf.open(filename, std::ios::binary);
    conf.read(reinterpret_cast<char *>(&v[0]), lsites * doublesPerSite * sizeof(double));
    conf.close();

    for (int site_index = 0; site_index < lsites; site_index++)
    {
        k = site_index * doublesPerSite;
        j = 0;
        for (int i = 0; i < 5; i++)
        {
            lattice[site_index].field[i](0, 0) = std::complex<double>(v[k + j++], v[k + j++]);
            lattice[site_index].field[i](0, 1) = std::complex<double>(v[k + j++], v[k + j++]);
            lattice[site_index].field[i](1, 0) = std::complex<double>(v[k + j++], v[k + j++]);
            lattice[site_index].field[i](1, 1) = std::complex<double>(v[k + j++], v[k + j++]);
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