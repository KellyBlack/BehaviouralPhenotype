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
    std::cout << "Starting" << std::endl;

    // Define the number of grid points to use.
    // Initialize it so that the matrices that are used to create the system
    // of equations are initialized.
    // Next, set the values of various constants.
    int N = 40;
    Butterflies theButterflies(N);
    theButterflies.setMu(1.0);
    theButterflies.setC(3.0);
    theButterflies.setF(2.0);
    theButterflies.setG(0.2);
    theButterflies.setD(0.5);

    // Variables used to save the results of calculations into a
    // data file.
    std::ofstream resultsFile;
    resultsFile.open(OUTPUTFILE);


    std::cout << "Pre-processing" << std::endl;
    // Define the Gauss quadrature.
    theButterflies.initializeLegendreParams();
    theButterflies.writeAbscissa(resultsFile);

    // Build the system and solve.
    std::cout << "Calculating an approximation" << std::endl;
    theButterflies.copyCurrentStateToTemp();
    bool canInvert(true);
    do
    {
        theButterflies.buildJacobian();
        canInvert = theButterflies.solveLinearizedSystem();
        if(canInvert){
            theButterflies.updateNewtonStep();
            std::cout << "  stepping" << std::endl;
        }

    } while((theButterflies.normDelta()>0.001) && canInvert);
    theButterflies.writeCurrentApprox(resultsFile);

    // Clean up the data file and close it
    resultsFile.close();

    std::cout << "Done" << std::endl;
    return(0);
}
