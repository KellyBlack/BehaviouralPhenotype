#include <iostream>
#include <fstream>
#include <cmath>

#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"
#include "butterflies.h"

#define OUTPUTFILE "approximation.csv"

int main()
{
    std::cout << "Pre-processing starting" << std::endl;

    // Define the number of grid points to use.
    // Define the matrices that are used to create the system of equations.
    int N = 40;
    butterflies theButterflies(N);

    // Variables used to save the results of calculations into a
    // data file.
    std::ofstream resultsFile;
    resultsFile.open(OUTPUTFILE);


    std::cout << "Pre-processing" << std::endl;
    // Define the Gauss quadrature.
    theButterflies.initializeLegendreParams();
    theButterflies.writeAbscissa(resultsFile);

    // Build the system.
    theButterflies.buildJacobian();
    theButterflies.solveLinearizedSystem();
    theButterflies.writeCurrentApprox(resultsFile);

    // Clean up the data file and close it
    resultsFile.close();

    std::cout << "Done" << std::endl;
    return(0);
}
