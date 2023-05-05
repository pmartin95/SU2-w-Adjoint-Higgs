#include "test.hpp"
#include <iomanip>
#include "hmc.hpp"
#include "lattice_ops.hpp"
void simulation1(int argc, char **argv)
{
    // arguments
    bool thermalize = true;
    if (argc > 0)
    {
        for (int i = 1; i < argc; i++)
        {
            if (!strcmp(argv[i], "-beta"))
            {
                beta = atof(argv[i + 1]);
            }
            if (!strcmp(argv[i], "-kappa"))
            {
                kappa = atof(argv[i + 1]);
            }
            if (!strcmp(argv[i], "-conf"))
            {
                thermalize = false;
                std::cout << "Skipping thermalization.\n";
                pullConfig("configurations/conf" + std::to_string(beta) + ".bin");
            }
        }
    }
    // File setup
    int highest_len = 3;
    // Setup data holding vectors
    std::vector<std::vector<double>> MxMdata, MxNdata;
    for (int i = 0; i < highest_len; i++) // Init data holding vectors
    {
        std::vector<double> tempdata1, tempdata2;
        MxMdata.push_back(std::move(tempdata1));
        MxNdata.push_back(std::move(tempdata2));
    }
    std::vector<double> MxMave(highest_len, 0.0), MxMerror(highest_len, 0.0), MxNave(highest_len, 0.0), MxNerror(highest_len, 0.0);

    // Thermalization
    if (thermalize)
    {
        hotLattice();
        for (int i = 0; i < 5000; i++)
        {
            metropolisHastingsSweep();
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
            }
            for (int i = 0; i < highest_len; i++)
            {
                MxMdata[i].push_back(rectangleAverage(i + 1, i + 1));
                MxNdata[i].push_back(rectangleAverage(i + 1, i));
            }
        }
        // Write configuration to file
        pushConfig("configurations/conf" + std::to_string(beta) + ".bin");
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
        // std::cout << "beta: " << beta << " X(1,1) = " << MxMave[0] << "\u00b1" << MxMerror[0] << " X(2,2)%=" << std::abs(MxMerror[1] / MxMave[1]) << " X(1,2)%=" << std::abs(MxNerror[1] / MxNave[1]) << std::endl;
    } while (check);

    // Rectangle file handling---------
    std::vector<std::string> MxMfilenames, MxNfilenames;
    std::vector<std::ofstream> MxMfiles, MxNfiles;
    for (int i = 0; i < highest_len; i++) // Name files and open them
    {
        MxMfilenames.push_back(datFolder + "rect" + std::to_string(i + 1) + "x" + std::to_string(i + 1) + "beta" + std::to_string(beta) + ".txt");
        MxNfilenames.push_back(datFolder + "rect" + std::to_string(i + 1) + "x" + std::to_string(i) + "beta" + std::to_string(beta) + ".txt");
        std::ofstream temp1, temp2;
        temp1.open(MxMfilenames[i], std::ios_base::app);
        temp2.open(MxNfilenames[i], std::ios_base::app);
        MxMfiles.push_back(std::move(temp1));
        MxNfiles.push_back(std::move(temp2));
    }

    for (int i = 0; i < highest_len; i++)
    {
        MxMfiles[i] << beta << " " << MxMave[i] << " " << MxMerror[i] << "\n";
        MxNfiles[i] << beta << " " << MxNave[i] << " " << MxNerror[i] << "\n";
    }
    for (auto &file : MxMfiles)
        file.close();
    for (auto &file : MxNfiles)
        file.close();
    //--------------------------------
}

void confIdentical()
{
    hotLattice();
    site *templattice = new site[lsites];
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
void simulation2(int argc, char **argv)
{
    int numThermalizations = 1000, iterBetweenObs = 100, numObs = 100;

    std::vector<double> higgsSquareData, plaquetteData;
    std::vector<matrix> higgsData;

    // kappa = 0, lambda = 0, m^2 = 0, beta = 2.3
    kappa = 0.0;
    lambda = 0.0;
    m2 = 0.0;
    hotLattice();
    std::cout << "\\kappa=" << kappa << " ,lambda=" << lambda << " ,m^2=" << m2 << " ,beta=" << beta << std::endl;
    for (int i = 0; i < numThermalizations; i++)
        metropolisHastingsSweep();
    for (int i = 0; i < numObs; i++)
    {
        for (int j = 0; j < iterBetweenObs; j++)
        {
            metropolisHastingsSweep();
        }
        higgsSquareData.push_back(higgsSquareAverage());
        plaquetteData.push_back(plaquetteAverage());
        higgsData.push_back(higgsAverage());
    }
    std::cout << "<Tr\\phi^2>=" << average(higgsSquareData) << " <Tr\\phi>=" << average(higgsData) << " Plaquette=" << average(plaquetteData) << std::endl;
    higgsSquareData.clear();
    plaquetteData.clear();
    higgsData.clear();
    // kappa = 1, lambda = 0, m^2 = 0, beta = 2.3
    kappa = 1.0;
    hotLattice();
    std::cout << "\\kappa=" << kappa << " ,lambda=" << lambda << " ,m^2=" << m2 << " ,beta=" << beta << std::endl;
    for (int i = 0; i < numThermalizations; i++)
        metropolisHastingsSweep();
    for (int i = 0; i < numObs; i++)
    {
        for (int j = 0; j < iterBetweenObs; j++)
        {
            metropolisHastingsSweep();
        }
        higgsSquareData.push_back(higgsSquareAverage());
        plaquetteData.push_back(plaquetteAverage());
        higgsData.push_back(higgsAverage());
    }
    std::cout << "<Tr\\phi^2>=" << average(higgsSquareData) << " <Tr\\phi>=" << average(higgsData) << " Plaquette=" << average(plaquetteData) << std::endl;
    higgsSquareData.clear();
    plaquetteData.clear();
    higgsData.clear();
    // kappa = 0, lambda = 0.1, m^2 = -0.2, beta = 2.3
    kappa = 0.0;
    lambda = 0.1;
    m2 = -0.2;
    hotLattice();
    std::cout << "\\kappa=" << kappa << " ,lambda=" << lambda << " ,m^2=" << m2 << " ,beta=" << beta << std::endl;
    for (int i = 0; i < numThermalizations; i++)
        metropolisHastingsSweep();
    for (int i = 0; i < numObs; i++)
    {
        for (int j = 0; j < iterBetweenObs; j++)
        {
            metropolisHastingsSweep();
        }
        higgsSquareData.push_back(higgsSquareAverage());
        plaquetteData.push_back(plaquetteAverage());
        higgsData.push_back(higgsAverage());
    }
    std::cout << "<Tr\\phi^2>=" << average(higgsSquareData) << " <Tr\\phi>=" << average(higgsData) << " Plaquette=" << average(plaquetteData) << std::endl;
    higgsSquareData.clear();
    plaquetteData.clear();
    higgsData.clear();
    // kappa = 1, lambda = 0.1, m^2 = -0.2, beta = 2.3
    kappa = 1.0;
    hotLattice();
    std::cout << "\\kappa=" << kappa << " ,lambda=" << lambda << " ,m^2=" << m2 << " ,beta=" << beta << std::endl;
    for (int i = 0; i < numThermalizations; i++)
        metropolisHastingsSweep();
    for (int i = 0; i < numObs; i++)
    {
        for (int j = 0; j < iterBetweenObs; j++)
        {
            metropolisHastingsSweep();
        }
        higgsSquareData.push_back(higgsSquareAverage());
        plaquetteData.push_back(plaquetteAverage());
        higgsData.push_back(higgsAverage());
    }
    std::cout << "<Tr\\phi^2>=" << average(higgsSquareData) << " <Tr\\phi>=" << average(higgsData) << " Plaquette=" << average(plaquetteData) << std::endl;
    higgsSquareData.clear();
    plaquetteData.clear();
    higgsData.clear();
}

// Observables to collect:
// correlator,higgsSquareAverage,
void simulation3(int argc, char **argv)
{
    std::vector<std::vector<double>> correlationData;
    for (int i = 0; i < ldir[0]; i++)
    {
        std::vector<double> temp;
        correlationData.push_back(std::move(temp));
    }

    std::vector<double> higgsSquareData;

    bool thermalize = true;
    if (argc > 0)
    {
        for (int i = 1; i < argc; i++)
        {
            if (!strcmp(argv[i], "-beta"))
            {
                beta = atof(argv[i + 1]);
                std::cout << "Using beta" << beta << std::endl;
            }
            if (!strcmp(argv[i], "-m2"))
            {
                m2 = atof(argv[i + 1]);
            }
            if (!strcmp(argv[i], "-conf"))
            {
                thermalize = false;
                std::cout << "Skipping thermalization.\n";
                pullConfig(confFolder + std::to_string(m2) + ".bin");
            }
        }
    }
    int numThermalizations = 3000, numObservations = 100, numSweepsPerObservation = 1000;
    if (thermalize)
    {
        hotLattice();
        for (int i = 0; i < numThermalizations; i++)
        {
            metropolisHastingsSweep();
        }
        std::cout << "Thermalized.\n";
    }

    for (int i = 0; i < numObservations; i++)
    {
        for (int j = 0; j < numSweepsPerObservation; j++)
        {
            metropolisHastingsSweep();
        }
        // collect observations
        for (int i = 0; i < ldir[0]; i++)
        {
            correlationData[i].push_back(averageCorrelatorVolume(i));
        }
        higgsSquareData.push_back(higgsSquareAverage());
        // write configuration to file
        pushConfig(confFolder + std::to_string(m2) + ".bin");
    }
    // calculate statistics
    std::vector<double> correlatorAggregate, correlatorError;
    double higgsSquareAggregate, higgsSquareError;

    // for (int i = 0; i < ldir[0]; i++)
    // {
    //     double tempAve, tempError;
    //     computeJackknifeStatistics(correlationData[i], 10, tempAve, tempError);
    //     correlatorAggregate.push_back(tempAve);
    //     correlatorError.push_back(tempError);
    // }
    // computeJackknifeStatistics(higgsSquareData, 10, higgsSquareAggregate, higgsSquareError);
    //! Commenting out this section temporarily
    //  write observations to file
    std::ofstream corrFile(datFolder + "correlation-m2-" + std::to_string(m2) + ".txt", std::ios_base::app);
    std::ofstream higgsFile(datFolder + "higgssquare-m2-" + std::to_string(m2) + ".txt", std::ios_base::app);
    for (int i = 0; i < ldir[0]; i++)
    {
        for (int j = 0; j < correlationData[i].size(); j++)
        {
            corrFile << m2 << " " << i << " " << correlationData[i][j] << std::endl;
        }
    }
    higgsFile << m2 << " " << higgsSquareAggregate << " " << higgsSquareError << std::endl;
}

// simulation 1 is supposed to collect rectangle data

void simulation1mod(int argc, char **argv)
{
    // File setup
    int highest_len = 3;

    // Setup data holding vectors
    std::vector<std::vector<double>> MxMdata, MxNdata;

    MxMdata.resize(highest_len, std::vector<double>(0));
    MxNdata.resize(highest_len, std::vector<double>(0));
    std::vector<double> MxMave(highest_len, 0.0), MxMerror(highest_len, 0.0), MxNave(highest_len, 0.0), MxNerror(highest_len, 0.0);

    // Thermalization
    if (thermalize)
    {
        hotLattice();
        for (int i = 0; i < 3000; i++)
        {
            metropolisHastingsSweep();
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
            polyakovLines("polyakov" + identifier + ".txt", 1, 2);
        }

        // Write configuration to file
        pushConfig("configurations/conf" + identifier + ".bin");
        // Statistics check for convergence
        for (int i = 0; i < highest_len; i++)
        {
            computeJackknifeStatistics(MxMdata[i], 10, MxMave[i], MxMerror[i]);
            computeJackknifeStatistics(MxNdata[i], 10, MxNave[i], MxNerror[i]);
        }
        // for (auto &dat : MxMdata[0])
        // {
        //     std::cout << dat << std::endl;
        // }
        // check = (std::abs(MxMerror[0] / MxMave[0]) > 0.05 ||
        //          std::abs(MxMerror[1] / MxMave[1]) > 0.1 ||
        //          std::abs(MxNerror[1] / MxNave[1]) > 0.1) &&
        //         (iter_count < MAX_ITER);
        check = (std::abs(MxMerror[0] / MxMave[0]) > 0.05) &&
                (iter_count < MAX_ITER);
    } while (check);

    // Rectangle file handling---------
    std::vector<std::string> MxMfilenames, MxNfilenames;
    std::vector<std::ofstream> MxMfiles, MxNfiles;
    std::cout << "num obs: " << MxMdata[0].size() << std::endl;
    for (int i = 0; i < highest_len; i++) // Name files and open them
    {
        MxMfilenames.push_back(datFolder + "rect" + std::to_string(i + 1) + "x" + std::to_string(i + 1) + identifier + ".txt");
        MxNfilenames.push_back(datFolder + "rect" + std::to_string(i + 1) + "x" + std::to_string(i) + identifier + ".txt");
        std::ofstream temp1, temp2;
        temp1.open(MxMfilenames[i], std::ios_base::app);
        temp2.open(MxNfilenames[i], std::ios_base::app);
        MxMfiles.push_back(std::move(temp1));
        MxNfiles.push_back(std::move(temp2));
    }

    for (int i = 0; i < highest_len; i++)
    {
        MxMfiles[i] << beta << " " << MxMave[i] << " " << MxMerror[i] << "\n";
        MxNfiles[i] << beta << " " << MxNave[i] << " " << MxNerror[i] << "\n";
    }
    for (auto &file : MxMfiles)
        file.close();
    for (auto &file : MxNfiles)
        file.close();
    //--------------------------------
}

void simulation2mod(int argc, char **argv)
{
    int numThermalizations = 1000, iterBetweenObs = 100, numObs = 100;

    std::vector<double> higgsSquareData, plaquetteData;
    std::vector<matrix> higgsData;

    // kappa = 0, lambda = 0, m^2 = 0, beta = 2.3
    kappa = 0.0;
    lambda = 0.0;
    m2 = 0.0;
    hotLattice();
    std::cout << "\\kappa=" << kappa << " ,lambda=" << lambda << " ,m^2=" << m2 << " ,beta=" << beta << std::endl;
    for (int i = 0; i < numThermalizations; i++)
        metropolisHastingsSweep();
    for (int i = 0; i < numObs; i++)
    {
        for (int j = 0; j < iterBetweenObs; j++)
        {
            metropolisHastingsSweep();
        }
        higgsSquareData.push_back(higgsSquareAverage());
        plaquetteData.push_back(plaquetteAverage());
        higgsData.push_back(higgsAverage());
    }
    std::cout << "<Tr\\phi^2>=" << average(higgsSquareData) << " <Tr\\phi>=" << average(higgsData) << " Plaquette=" << average(plaquetteData) << std::endl;
    higgsSquareData.clear();
    plaquetteData.clear();
    higgsData.clear();
    // kappa = 1, lambda = 0, m^2 = 0, beta = 2.3
    kappa = 1.0;
    hotLattice();
    std::cout << "\\kappa=" << kappa << " ,lambda=" << lambda << " ,m^2=" << m2 << " ,beta=" << beta << std::endl;
    for (int i = 0; i < numThermalizations; i++)
        metropolisHastingsSweep();
    for (int i = 0; i < numObs; i++)
    {
        for (int j = 0; j < iterBetweenObs; j++)
        {
            metropolisHastingsSweep();
        }
        higgsSquareData.push_back(higgsSquareAverage());
        plaquetteData.push_back(plaquetteAverage());
        higgsData.push_back(higgsAverage());
    }
    std::cout << "<Tr\\phi^2>=" << average(higgsSquareData) << " <Tr\\phi>=" << average(higgsData) << " Plaquette=" << average(plaquetteData) << std::endl;
    higgsSquareData.clear();
    plaquetteData.clear();
    higgsData.clear();
    // kappa = 0, lambda = 0.1, m^2 = -0.2, beta = 2.3
    kappa = 0.0;
    lambda = 0.1;
    m2 = -0.2;
    hotLattice();
    std::cout << "\\kappa=" << kappa << " ,lambda=" << lambda << " ,m^2=" << m2 << " ,beta=" << beta << std::endl;
    for (int i = 0; i < numThermalizations; i++)
        metropolisHastingsSweep();
    for (int i = 0; i < numObs; i++)
    {
        for (int j = 0; j < iterBetweenObs; j++)
        {
            metropolisHastingsSweep();
        }
        higgsSquareData.push_back(higgsSquareAverage());
        plaquetteData.push_back(plaquetteAverage());
        higgsData.push_back(higgsAverage());
    }
    std::cout << "<Tr\\phi^2>=" << average(higgsSquareData) << " <Tr\\phi>=" << average(higgsData) << " Plaquette=" << average(plaquetteData) << std::endl;
    higgsSquareData.clear();
    plaquetteData.clear();
    higgsData.clear();
    // kappa = 1, lambda = 0.1, m^2 = -0.2, beta = 2.3
    kappa = 1.0;
    hotLattice();
    std::cout << "\\kappa=" << kappa << " ,lambda=" << lambda << " ,m^2=" << m2 << " ,beta=" << beta << std::endl;
    for (int i = 0; i < numThermalizations; i++)
        metropolisHastingsSweep();
    for (int i = 0; i < numObs; i++)
    {
        for (int j = 0; j < iterBetweenObs; j++)
        {
            metropolisHastingsSweep();
        }
        higgsSquareData.push_back(higgsSquareAverage());
        plaquetteData.push_back(plaquetteAverage());
        higgsData.push_back(higgsAverage());
    }
    std::cout << "<Tr\\phi^2>=" << average(higgsSquareData) << " <Tr\\phi>=" << average(higgsData) << " Plaquette=" << average(plaquetteData) << std::endl;
    higgsSquareData.clear();
    plaquetteData.clear();
    higgsData.clear();
}

// Observables to collect:
// correlator,higgsSquareAverage,
void simulation3mod(int argc, char **argv)
{
    std::vector<std::vector<double>> correlationData;
    for (int i = 0; i < ldir[0]; i++)
    {
        std::vector<double> temp;
        correlationData.push_back(std::move(temp));
    }

    std::vector<double> higgsSquareData;
    int numThermalizations = configureStep(), numObservations = 100, numSweepsPerObservation = numThermalizations / 3;

    if (thermalize)
    {
        hotLattice();
        for (int i = 0; i < numThermalizations; i++)
        {
            metropolisHastingsSweep();
        }
        std::cout << "Thermalized.\n";
    }

    for (int i = 0; i < numObservations; i++)
    {
        for (int j = 0; j < numSweepsPerObservation; j++)
        {
            metropolisHastingsSweep();
        }
        // collect observations
        for (int i = 0; i < ldir[0]; i++)
        {
            correlationData[i].push_back(averageCorrelatorVolume(i));
        }
        higgsSquareData.push_back(higgsSquareAverage());
        // write configuration to file
        pushConfig(confFolder + std::to_string(m2) + ".bin");
    }
    // calculate statistics
    std::vector<double> correlatorAggregate, correlatorError;
    double higgsSquareAggregate, higgsSquareError;

    // for (int i = 0; i < ldir[0]; i++)
    // {
    //     double tempAve, tempError;
    //     computeJackknifeStatistics(correlationData[i], 10, tempAve, tempError);
    //     correlatorAggregate.push_back(tempAve);
    //     correlatorError.push_back(tempError);
    // }
    // computeJackknifeStatistics(higgsSquareData, 10, higgsSquareAggregate, higgsSquareError);
    //! Commenting out this section temporarily
    //  write observations to file
    std::ofstream corrFile(datFolder + "correlation-m2-" + std::to_string(m2) + ".txt", std::ios_base::app);
    std::ofstream higgsFile(datFolder + "higgssquare-m2-" + std::to_string(m2) + ".txt", std::ios_base::app);
    for (int i = 0; i < ldir[0]; i++)
    {
        for (int j = 0; j < static_cast<int>(correlationData[i].size()); j++)
        {
            corrFile << m2 << " " << i << " " << correlationData[i][j] << std::endl;
        }
    }
    higgsFile << m2 << " " << higgsSquareAggregate << " " << higgsSquareError << std::endl;
    corrFile.close();
    higgsFile.close();
}

void simulation4(int argc, char **argv)
{
    //--------Data to declare
    //====================================================================================================
    // cold/hot plaquette/Higgs square thermalization convergence
    // rectangles
    // Higgs square
    // polyakov loops
    // correlation functions

    // Thermalization
    //====================================================================================================
    int nTherm = 3000;
    int simKey = 1;
    while (fileExists(datFolder + "higgsSquare" + identifier + std::to_string(simKey) + ".dat"))
    {
        simKey++;
    }
    if (thermalize)
    {
        std::ofstream hotThermPlaq(datFolder + "hotThermPlaq" + identifier + std::to_string(simKey) + ".dat", std::ios_base::app);
        std::ofstream hotThermHiggsSquare(datFolder + "hotThermHiggsSquare" + identifier + std::to_string(simKey) + ".dat", std::ios_base::app);
        std::ofstream coldThermPlaq(datFolder + "coldThermPlaq" + identifier + std::to_string(simKey) + ".dat", std::ios_base::app);
        std::ofstream coldThermHiggsSquare(datFolder + "coldThermHiggsSquare" + identifier + std::to_string(simKey) + ".dat", std::ios_base::app);
        // nTherm = configureStep();
        nTherm = 3000;
        coldLattice();
        for (int i = 0; i < nTherm; i++)
        {
            metropolisHastingsSweep();
            coldThermPlaq << plaquetteAverage() << std::endl;
            coldThermHiggsSquare << higgsSquareAverage() << std::endl;
        }
        coldThermPlaq.close();
        coldThermHiggsSquare.close();
        hotLattice();
        for (int i = 0; i < nTherm; i++)
        {
            metropolisHastingsSweep();
            hotThermPlaq << plaquetteAverage() << std::endl;
            hotThermHiggsSquare << higgsSquareAverage() << std::endl;
        }
        hotThermPlaq.close();
        hotThermHiggsSquare.close();
    }
    else if (configSteps)
    {
        // nTherm = configureStep();
        nTherm = 3000;
    }
    std::ofstream sim_info(datFolder + "simInfo" + identifier + ".txt");
    sim_info << "nTherm= " << nTherm << std::endl;
    sim_info << "Higgs Step= " << step_size_higgs << std::endl;
    sim_info << "Link Step= " << rot_size << std::endl;
    std::cout << "Hello.\n";
    sim_info.close();

    // Sweeps
    //====================================================================================================
    int nIter = nTherm / 3;
    int nObs = 1000;
    std::cout << "Sweeping.\n";
    for (int i = 0; i < nObs; i++)
    {
        for (int j = 0; j < nIter; j++)
        {
            metropolisHastingsSweep();
        }
        // Write config to file
        pushConfig(confFolder + identifier + ".bin");
        std::ofstream data_out;
        //  Collect data
        //  rectangles
        data_out.open(datFolder + "wilsonRectangles" + identifier + std::to_string(simKey) + ".dat", std::ios_base::app);
        for (int i = 1; i < 4; i++)
        {
            for (int j = 1; j <= i; j++)
            {
                data_out << beta << " " << m2 << " " << i << " " << j << " " << rectangleAverage(i, j) << std::endl;
            }
        }
        data_out.close();
        //  Higgs square average
        data_out.open(datFolder + "higgsSquare" + identifier + std::to_string(simKey) + ".dat", std::ios_base::app);
        data_out << beta << " " << m2 << " " << higgsSquareAverage() << std::endl;
        data_out.close();

        // Higgs Square Sum
        data_out.open(datFolder + "higgsSquareSum" + identifier + std::to_string(simKey) + ".dat", std::ios_base::app);
        data_out << beta << " " << m2 << " " << higgsSquareSum() << std::endl;
        data_out.close();
        //  polyakov loops
        polyakovLines("polyakovLoops" + identifier + std::to_string(simKey) + ".dat", 1, 2);
        //  correlation functions
        data_out.open(datFolder + "higgsCorrelator" + identifier + std::to_string(simKey) + ".dat", std::ios_base::app);
        for (int t = 0; t < lt; t++)
        {
            data_out << beta << " " << m2 << " " << t << " " << averageCorrelatorVolume(t) << std::endl;
        }
        data_out.close();
    }
}

void argumentInput(int argc, char **argv)
{
    if (argc > 0)
    {
        for (int i = 1; i < argc; i++)
        {
            if (!strcmp(argv[i], "-beta"))
            {
                beta = atof(argv[i + 1]);
                std::cout << "Using beta " << beta << std::endl;
            }
            else if (!strcmp(argv[i], "-m2"))
            {
                m2 = atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-lambda"))
            {
                lambda = atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-kappa"))
            {
                kappa = atof(argv[i + 1]);
            }
            else if (!strcmp(argv[i], "-llength"))
            {
                l = atoi(argv[i + 1]);
                ldir[1] = l;
                ldir[2] = l;
                ldir[3] = l;
            }
            else if (!strcmp(argv[i], "-ltime"))
            {
                lt = atoi(argv[i + 1]);
                ldir[0] = lt;
            }
            else if (!strcmp(argv[i], "-bc"))
            {
                if (!strcmp(argv[i + 1], "p"))
                {
                    bc = &ptwist;
                    bcName = 'p';
                }
                else if (!strcmp(argv[i + 1], "c"))
                {
                    bc = &ctwist;
                    bcName = 'c';
                }
                else if (!strcmp(argv[i + 1], "t"))
                {
                    bc = &twist;
                    bcName = 't';
                }
            }
        }
        lsites = lt * l * l * l;
        createIdentifier(identifier);
        for (int i = 1; i < argc; i++)
        {
            if (!strcmp(argv[i], "-conf"))
            {
                thermalize = false;
                std::cout << "Using previous configuration instead of thermalization.\n";
                pullConfig(confFolder + identifier + ".bin");
            }
        }
    }
    lattice = new site[lsites];
    lattice1_global = new site[lsites];
    lattice2_global = new site[lsites];
    plattice = new site[lsites];
    plattice1_global = new site[lsites];
    plattice2_global = new site[lsites];

    std::random_device r;
    rng.seed(r());
}

int configureStep()
{
    int higgsDistance = 250, linkDistance = 1400;
    bool statusH = true, statusL = true;
    double higgsAcceptance, linkAcceptance;
    int maxI = 11, iter = 0;
    do
    {
        iter++;
        Naccept = 0;
        Nreject = 0;
        NacceptLink = 0;
        NrejectLink = 0;
        NacceptHiggs = 0;
        NrejectHiggs = 0;

        hotLattice();
        for (int i = 0; i < 1000; i++)
        {
            metropolisHastingsSweep();
        }
        higgsAcceptance = static_cast<double>(NacceptHiggs) / 1000;
        linkAcceptance = static_cast<double>(NacceptLink) / 1000;

        if (higgsAcceptance > 0.8)
        {
            step_size_higgs += 0.03;
        }
        else if (higgsAcceptance < 0.6)
        {
            step_size_higgs -= 0.03;
        }
        else
        {
            statusH = false;
        }
        if (linkAcceptance > 0.8)
        {
            rot_size += 0.03;
        }
        else if (linkAcceptance < 0.6)
        {
            rot_size -= 0.03;
        }
        else
        {
            statusL = false;
        }

    } while ((statusH || statusL) && (iter < maxI));
    double fIter = static_cast<double>(14) / 0.7 / rot_size;
    return 150 * static_cast<int>(ceil(fIter));
}

void boundaryConditionTest(int argc, char **argv)
{
    // p bc tests
    bc = &ptwist;
    hotLattice();
    for (int site_index = 0; site_index < lsites; site_index++)
        for (int i = 0; i < 4; i++)
        {
            int jump[4] = {0};
            jump[i] = ldir[i];
            for (int j = 0; j < 5; j++)
            {
                // if (!verify::equivalent(callLatticeSite(site_index, jump, j), callLatticeSite(site_index, jumpNone, j)))
                //     std::cout << "problem.\n";
                if (!verify::equivalent(callLatticeSite(site_index, jump, j), lattice[site_index].field[j]))
                {
                    // std::cout << "problem.\n";
                    // std::cout << callLatticeSite(site_index, jump, j) << "\t" << lattice[site_index].field[j] << std::endl;
                    std::cout << site_index << std::endl;
                }
            }
        }

    // c twist tests
    // bc = &ctwist;
    // hotLattice();
    // for (int site_index = 0; site_index < lsites; site_index++)
    //     for (int i = 0; i < 4; i++)
    //     {
    //         int jump[4] = {0};
    //         int jumpNone[4] = {};
    //         jump[i] = ldir[i];
    //         for (int j = 0; j < 5; j++)
    //         {
    //             if (!verify::equivalent(callLatticeSite(site_index, jump, j), callLatticeSite(site_index, jumpNone, j)))
    //                 std::cout << "problem.\n";
    //             matrix temp = lattice[site_index].field[j];
    //             ctwist(temp, j, i, 1);
    //             if (!verify::equivalent(callLatticeSite(site_index, jump, j), temp))
    //                 std::cout << "problem.\n";
    //         }
    //     }
}

void indexTest()
{
    for (int i = 0; i < lsites; i++)
    {
        int r[4], j;
        siteIndexToCoordinates(i, r);
        j = coordinatesToSiteIndex(r);
        if (i != j)
            std::cout << i << " " << j << std::endl;
    }
}

void simulation5(int argc, char **argv)
{
    // Measure Polaykov loops, Wilson Rectangles, Higgs fields, Higgs correlators
    // Thermalization
    std::vector<std::vector<double>> MxMdata, MxNdata;
    int highest_len = 3; // Highest length of Polyakov loop that is checked against for convergence
    std::vector<double> MxMave(highest_len, 0.0), MxMerror(highest_len, 0.0), MxNave(highest_len, 0.0), MxNerror(highest_len, 0.0);

    for (int i = 0; i < highest_len; i++) // Init data holding vectors
    {
        std::vector<double> tempdata1, tempdata2;
        MxMdata.push_back(std::move(tempdata1));
        MxNdata.push_back(std::move(tempdata2));
    }
    if (thermalize)
    {
        hotLattice();
        for (int i = 0; i < 3000; i++)
        {
            metropolisHastingsSweep();
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
            }
            midSimObservables();
            for (int i = 0; i < highest_len; i++)
            {
                MxMdata[i].push_back(rectangleAverage(i + 1, i + 1));
                MxNdata[i].push_back(rectangleAverage(i + 1, i));
            }
        }
        // Write configuration to file
        pushConfig(confFolder + identifier + ".bin");
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
    } while (check);
}

void midSimObservables()
{
    // Measure Polaykov loops, Wilson Rectangles, Higgs fields, Higgs correlators
    polyakovLines(datFolder + "polyakov" + identifier + ".dat", 1, 2);
    for (int t = 1; t < lt; t++)
    {
        std::ofstream file(datFolder + "rect" + std::to_string(t) + "x" + std::to_string(t) + identifier + ".dat", std::ios::app);
        file << std::setprecision(16) << std::fixed << std::scientific;
        file << rectangleAverage(t, t) << std::endl;
    }
    for (int t = 2; t < lt; t++)
    {
        std::ofstream file(datFolder + "rect" + std::to_string(t) + "x" + std::to_string(t - 1) + identifier + ".dat", std::ios::app);
        file << std::setprecision(16) << std::fixed << std::scientific;
        file << rectangleAverage(t, t - 1) << std::endl;
    }
    std::ofstream higgsFile(datFolder + "higgsSquare" + identifier + ".dat", std::ios::app);
    higgsFile << std::setprecision(16) << std::fixed << std::scientific;
    higgsFile << higgsSquareAverage() << std::endl;
    higgsFile.close();

    std::ofstream higgsCorrFile(datFolder + "higgsCorr" + identifier + ".dat", std::ios::app);
    for (int t = 0; t < lt; t++)
    {
        higgsCorrFile << t << " ";
        higgsCorrFile << std::setprecision(16) << std::fixed << std::scientific;
        higgsCorrFile << averageCorrelatorVolume(t) << "\n";
    }
    higgsCorrFile.close();
}

void simulationHMC1(int argc, char **argv)
{
    // Initialize variables
    int sweep_between_obs = 40;
    int obs_between_checks = 100;
    int highest_len = 3; // Highest length of Polyakov loop that is checked against for convergence
    int iter_count = 0;
    bool check;

    // Measure Polaykov loops, Wilson Rectangles, Higgs fields, Higgs correlators
    // Initialize data structures for analysis
    std::vector<std::vector<double>> MxMdata, MxNdata;
    std::vector<double> MxMave(highest_len, 0.0), MxMerror(highest_len, 0.0), MxNave(highest_len, 0.0), MxNerror(highest_len, 0.0);

    for (int i = 0; i < highest_len; i++) // Init data holding vectors
    {
        std::vector<double> tempdata1, tempdata2;
        MxMdata.push_back(std::move(tempdata1));
        MxNdata.push_back(std::move(tempdata2));
    }
    hotLattice();
    // Thermalization
    if (thermalize)
    {
        hotLattice();
        for (int i = 0; i < 50; i++)
        {
            HMC_warmup(20);
        }
    }
    // Thermalization

    // Main simulation loop
    do
    {
        ++iter_count;

        // Perform sweeps and collect observables
        for (int i = 0; i < obs_between_checks; i++)
        {
            for (int j = 0; j < sweep_between_obs; j++)
            {
                // std::cout << WilsonAction(lattice) << std::endl;
                HMC(1000);
                std::cout << plaquetteAverage() << std::endl;
            }
            std::cout << "Acceptance rate: " << static_cast<double>(Naccept) / static_cast<double>(Naccept + Nreject) << std::endl;
            Naccept = 0;
            Nreject = 0;
            // midSimObservables(); //! pasting this out for now. All I want to do is keep track of plaquettes
            for (int i = 0; i < highest_len; i++)
            {
                MxMdata[i].push_back(rectangleAverage(i + 1, i + 1));
                MxNdata[i].push_back(rectangleAverage(i + 1, i));
            }
        }

        // Save configuration to file
        pushConfig(confFolder + identifier + ".bin");

        // Check for convergence
        for (int i = 0; i < highest_len; i++)
        {
            computeJackknifeStatistics(MxMdata[i], 10, MxMave[i], MxMerror[i]);
            computeJackknifeStatistics(MxNdata[i], 10, MxNave[i], MxNerror[i]);
        }

        check = (std::abs(MxMerror[0] / MxMave[0]) > 0.05 || std::abs(MxMerror[1] / MxMave[1]) > 0.1 || std::abs(MxNerror[1] / MxNave[1]) > 0.1) && (iter_count < MAX_ITER);
    } while (check);
}

void testExpCK(int num_tests)
{
    std::mt19937 gen(42);                               // Seed the random number generator with a constant for reproducibility
    std::uniform_real_distribution<double> dist(-1, 1); // Uniform distribution between -1 and 1

    for (int test = 0; test < num_tests; test++)
    {
        // Generate a random traceless Hermitian matrix
        matrix H;
        double realDiagonal = dist(gen);
        H(0, 0) = std::complex<double>(realDiagonal, 0);
        H(1, 1) = std::complex<double>(-realDiagonal, 0);
        std::complex<double> offDiagonal(dist(gen), dist(gen));
        H(0, 1) = offDiagonal;
        H(1, 0) = std::conj(offDiagonal);

        // Generate a random constant

        double t = dist(gen);

        // Compute the exponentiated matrix
        matrix result;
        expCK(t * I * H, result);

        // Check if the output is SU(2)
        if (!isSU2(result))
        {
            std::cout << "Test failed for t = " << t << " and H = \n"
                      << H << std::endl;
            std::cout << "Result = \n"
                      << result << std::endl;
        }
        else
        {
            // std::cout << "Test passed for t = " << t << " and H = \n"
            //           << H << std::endl;
        }
    }
}