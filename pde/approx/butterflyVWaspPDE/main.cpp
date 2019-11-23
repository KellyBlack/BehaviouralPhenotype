#include <iostream>
#include <fstream>
#include <cmath>

#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"
#include "butterflies.h"

//#define OUTPUTFILE "approximation.csv"
#define BINARYOUTPUTFILE "approximation.bin"
#define SKIP_PRINT_UPDATE 8000
#define SKIP_FILE_SAVE 75
#define NUMBER_TIME_LOOP 1000000
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
    theButterflies.setMu(0.5);
    theButterflies.setC(0.7);
    theButterflies.setG(4.0);
    theButterflies.setD(2.0);
    theButterflies.setA(0.6 /*0.5*/ /*0.7*/);
    theButterflies.setDT(dt);


    // Variables used to save the results of calculations into a
    // data file.
#ifdef OUTPUTFILE
    std::ofstream resultsFile;
    resultsFile.open(OUTPUTFILE);
#endif
#ifdef BINARYOUTPUTFILE
    std::fstream binFile (BINARYOUTPUTFILE, std::ios::out | std::ios::binary);
#endif

    std::cout << "Pre-processing" << std::endl;
    // Define the Gauss quadrature.
    theButterflies.initializeLegendreParams();
    theButterflies.initializeButterflies();
#ifdef OUTPUTFILE
    theButterflies.writeAbscissa(resultsFile);
#endif
#ifdef BINARYOUTPUTFILE
    theButterflies.writeParameters(binFile);
    theButterflies.writeBinaryHeader(binFile);
#endif

    // Start the time loop, and calculation an approximation at
    // each time step.
    double stepDeltaNorm = 0.0;
    for(timeLupe=0;(timeLupe<NUMBER_TIME_LOOP)&&(stepDeltaNorm<MAX_DELTA_NORM);++timeLupe)
    {
        t = static_cast<double>(timeLupe)*dt;

        if(timeLupe%(SKIP_PRINT_UPDATE)==0)
        {
            std::cout << "Calculating an approximation: " << timeLupe << " (" << t << ")" << std::endl;
        }

        if(timeLupe%SKIP_FILE_SAVE==0)
        {
#ifdef OUTPUTFILE
            theButterflies.writeCurrentApprox(t,resultsFile);
#endif
#ifdef BINARYOUTPUTFILE
            theButterflies.writeBinaryCurrentApprox(t,binFile);
#endif
        }

        // Build the system and solve.
        theButterflies.calculateRHS();
        theButterflies.copyCurrentStateToTemp();
        bool canInvert(true);
        int maxNewtonSteps = MAX_NEWTON_STEPS;
        do
        {
            // Perform the Newton steps to approximate the nonlinear
            // equations associated with the implicit system.
            theButterflies.buildJacobian();
            canInvert = theButterflies.solveLinearizedSystem();
            if(canInvert)
            {
                theButterflies.updateNewtonStep();
                if(timeLupe%(SKIP_PRINT_UPDATE)==0)
                    std::cout << "  step: " << MAX_NEWTON_STEPS - maxNewtonSteps << ", ";
                stepDeltaNorm = theButterflies.normDelta();
            }
            else
            {
                std::cout << "  System not invertible" << std::endl;
                maxNewtonSteps = 0;
            }

        } while((stepDeltaNorm>MAX_DELTA_NORM) && canInvert && (maxNewtonSteps-- > 0));

    }
    if(timeLupe%(SKIP_PRINT_UPDATE)==0)
        std::cout << std::endl;


    // Clean up the data file and close it
#ifdef OUTPUTFILE
    resultsFile.close();
#endif
#ifdef BINARYOUTPUTFILE
    binFile.close();
#endif

    std::cout << "Done" << std::endl;
    return(0);
}
