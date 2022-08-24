#include "test.hpp"
void simulation1(int argc, char **argv)
{
    // arguments
    double input_beta = 2.3;
    bool thermalize = true;
    if (argc > 0)
    {
        for (int i = 1; i < argc; i++)
        {
            if (!strcmp(argv[i], "-beta"))
            {
                input_beta = atof(argv[i + 1]);
                std::cout << "Using beta" << input_beta << std::endl;
            }
            if (!strcmp(argv[i], "-conf"))
            {
                thermalize = false;
                std::cout << "Skipping thermalization.\n";
                pullConfig("configurations/conf" + std::to_string(beta) + ".bin");
            }
        }

        std::cout << "Running beta= " << input_beta << std::endl;
    }

    // File setup
    int highest_len = 3;
    std::string base_filename = "rect";
    std::vector<std::string> MxMfilenames, MxNfilenames;
    std::vector<std::ofstream> MxMfiles, MxNfiles;
    std::vector<double> rect2x2observe;
    for (int i = 0; i < highest_len; i++) // Name files and open them
    {
        MxMfilenames.push_back(base_filename + std::to_string(i + 1) + "x" + std::to_string(i + 1) + ".txt");
        MxNfilenames.push_back(base_filename + std::to_string(i + 1) + "x" + std::to_string(i) + ".txt");
        std::ofstream temp1, temp2;
        temp1.open(MxMfilenames[i], std::ios_base::app);
        temp2.open(MxNfilenames[i], std::ios_base::app);
        MxMfiles.push_back(std::move(temp1));
        MxNfiles.push_back(std::move(temp2));
    }

    // Setup Beta distribution 0.0 to 4.0 with 0.1 in between
    std::vector<double> beta_distribution; // Create beta distribution
    beta_distribution.push_back(input_beta);

    // Cycle through beta distribution
    for (auto &B : beta_distribution)
    {
        beta = B;
        // Setup data holding vectors
        std::vector<std::vector<double>> MxMdata, MxNdata;
        for (int i = 0; i < highest_len; i++) // Init data holding vectors
        {
            std::vector<double> tempdata1, tempdata2;
            MxMdata.push_back(std::move(tempdata1));
            MxNdata.push_back(std::move(tempdata2));
        }
        std::vector<double> MxMave(highest_len, 0.0), MxMerror(highest_len, 0.0), MxNave(highest_len, 0.0), MxNerror(highest_len, 0.0);
        double ave, error;

        // Thermalization
        if (thermalize)
        {
            hotLattice();
            for (int i = 0; i < 5000; i++)
            {
                metropolisHastingsSweep();
                polyakovLinesAbs("hotAbsPoly" + std::to_string(beta) + ".txt");
                // std::cout << rectangleAverage(2, 1) << std::endl;
                // std::cout << lattice[5].field[2] << std::endl;
                if (i % 100 == 0)
                    std::cout << "Acceptance rate in thermalization:" << static_cast<double>(Naccept) / static_cast<double>(Naccept + Nreject) << std::endl;
                Naccept = 0;
                Nreject = 0;
            }
        }

        // Sweeps
        iter_count = 0;
        int sweep_between_obs = 1000, obs_between_checks = 100;
        bool check;
        do
        {
            ++iter_count;
            for (int i = 0; i < obs_between_checks; i++)
            {
                for (int j = 0; j < sweep_between_obs; j++)
                {
                    metropolisHastingsSweep();
                    Naccept = 0;
                    Nreject = 0;
                }
                for (int i = 0; i < highest_len; i++)
                {
                    MxMdata[i].push_back(rectangleAverage(i + 1, i + 1));
                    MxNdata[i].push_back(rectangleAverage(i + 1, i));
                }
                polyakovLines("polyakov" + std::to_string(beta) + ".txt", 1, 2);
            }
            // Write configuration to file
            pushConfig("configurations/conf" + std::to_string(beta) + ".bin");
            std::cout << "Configuration written to file.\n";
            // Statistics check for convergence
            for (int i = 0; i < highest_len; i++)
            {
                computeJackknifeStatistics(MxMdata[i], 10, MxMave[i], MxMerror[i]);
                computeJackknifeStatistics(MxNdata[i], 10, MxNave[i], MxNerror[i]);
            }
            for (auto &dat : MxMdata[0])
            {
                std::cout << dat << std::endl;
            }
            check = (std::abs(MxMerror[0] / MxMave[0]) > 0.05 || std::abs(MxMerror[1] / MxMave[1]) > 0.1 || std::abs(MxNerror[1] / MxNave[1]) > 0.1) && (iter_count < MAX_ITER);
            std::cout << "beta: " << beta << " X(1,1) = " << MxMave[0] << "\u00b1" << MxMerror[0] << " X(2,2)%=" << std::abs(MxMerror[1] / MxMave[1]) << " X(1,2)%=" << std::abs(MxNerror[1] / MxNave[1]) << std::endl;
        } while (check);

        for (int i = 0; i < highest_len; i++)
        {
            MxMfiles[i] << beta << " " << MxMave[i] << " " << MxMerror[i] << "\n";
            MxNfiles[i] << beta << " " << MxNave[i] << " " << MxNerror[i] << "\n";
        }
    }

    delete[] lattice;
    for (auto &file : MxMfiles)
        file.close();
    for (auto &file : MxNfiles)
        file.close();
}

void confIdentical()
{
    hotLattice();
    link *templattice = new link[lsites];
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            templattice[site_index].field[dir] = lattice[site_index].field[dir];
        }
    }

    pushConfig("configurations/test_conf.bin");
    pullConfig("configurations/test_conf.bin");
    for (int site_index = 0; site_index < lsites; site_index++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    std::string word1, word2;
                    word1 = std::to_string(lattice[site_index].field[dir](i, j).real());
                    word2 = std::to_string(templattice[site_index].field[dir](i, j).real());
                    if (strcmp(word1.c_str(), word2.c_str()))
                        std::cout << "halt.\n";
                    word1 = std::to_string(lattice[site_index].field[dir](i, j).imag());
                    word2 = std::to_string(templattice[site_index].field[dir](i, j).imag());
                    if (strcmp(word1.c_str(), word2.c_str()))
                        std::cout << "halt.\n";
                }
            }
        }
    }

    delete[] templattice;
}