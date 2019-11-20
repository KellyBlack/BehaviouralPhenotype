#include <iostream>
#include <fstream>
#include <cmath>

#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"
#include "butterflies.h"

#define OUTPUTFILE "approximation.csv"
#define NUMBER_TIME_LOOP 1000
#define MAX_NEWTON_STEPS 50
#define LEGENDRE_POLY_DEGREE 30
#define MAX_DELTA_NORM 0.0001

int main()
{
    // Set up the temporal variables.
    double t        = 0.0;
    double dt       = 0.000185;
    int    timeLupe = 0;
    std::cout << "Starting" << std::endl;

    // Define the number of grid points to use.
    // Initialize it so that the matrices that are used to create the system
    // of equations are initialized.
    // Next, set the values of various constants.
    int N = LEGENDRE_POLY_DEGREE;
    Butterflies theButterflies(N,N+2);
    theButterflies.setMu(1.0);
    theButterflies.setC(3.0);
    theButterflies.setF(2.0);
    theButterflies.setG(0.2);
    theButterflies.setD(0.5);
    theButterflies.setDT(dt);

    // Variables used to save the results of calculations into a
    // data file.
    std::ofstream resultsFile;
    resultsFile.open(OUTPUTFILE);


    std::cout << "Pre-processing" << std::endl;
    // Define the Gauss quadrature.
    theButterflies.initializeLegendreParams();
    theButterflies.writeAbscissa(resultsFile);

    // Start the time loop, and calculation an approximation at
    // each time step.
    double stepDeltaNorm = 0.0;
    for(timeLupe=0;(timeLupe<NUMBER_TIME_LOOP)&&(stepDeltaNorm<MAX_DELTA_NORM);++timeLupe)
    {
        t = static_cast<double>(timeLupe)*dt;

        // Build the system and solve.
        std::cout << "Calculating an approximation" << std::endl;
        theButterflies.calculateRHS();
        theButterflies.copyCurrentStateToTemp();
        bool canInvert(true);
        int maxNewtonSteps = MAX_NEWTON_STEPS;
        do
        {
            theButterflies.buildJacobian();
            canInvert = theButterflies.solveLinearizedSystem();
            if(canInvert)
            {
                theButterflies.updateNewtonStep();
                std::cout << "  stepping " << MAX_NEWTON_STEPS - maxNewtonSteps << std::endl;
                stepDeltaNorm = theButterflies.normDelta();
            }
            else
            {
                std::cout << "  System not invertible" << std::endl;
            }

        } while((stepDeltaNorm>MAX_DELTA_NORM) && canInvert && (maxNewtonSteps-- > 0));
        theButterflies.writeCurrentApprox(t,resultsFile);

    }

    // Clean up the data file and close it
    resultsFile.close();

    std::cout << "Done" << std::endl;
    return(0);
}
