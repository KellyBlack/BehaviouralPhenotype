#include <iostream>
#include <fstream>
#include <cmath>

#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"
#include "butterflies.h"

//#define OUTPUTFILE "approximation.csv"
#define BINARYOUTPUTFILE "approximation.bin"
#define TIMESKIP 30
#define NUMBER_TIME_LOOP 100000000
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
    theButterflies.setMu(0.01); //(0.02); //(0.02); //(0.02); //0.1);
    theButterflies.setC(0.8); //(0.80); //(0.9); //(1.4); //0.9);
    //theButterflies.setF(); //2.0);
    theButterflies.setG(1.5); //(1.1); //(0.8); //(1.5); //1.8);
    theButterflies.setD(0.4); //(0.3); //(0.4); //(0.6); //0.6);
    theButterflies.setA(0.4); //(0.3); //(0.1); //(0.4);
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
#ifdef OUTPUTFILE
    theButterflies.writeAbscissa(resultsFile);
#endif
#ifdef BINARYOUTPUTFILE
    theButterflies.writeBinaryHeader(binFile);
#endif

    // Start the time loop, and calculation an approximation at
    // each time step.
    double stepDeltaNorm = 0.0;
    for(timeLupe=0;(timeLupe<NUMBER_TIME_LOOP)&&(stepDeltaNorm<MAX_DELTA_NORM);++timeLupe)
    {
        t = static_cast<double>(timeLupe)*dt;

        if(timeLupe%(TIMESKIP)==0)
        {
#ifdef OUTPUTFILE
            theButterflies.writeCurrentApprox(t,resultsFile);
#endif
#ifdef BINARYOUTPUTFILE
            theButterflies.writeBinaryCurrentApprox(t,binFile);
#endif
            std::cout << "Calculating an approximation: " << timeLupe << " (" << t << ")" << std::endl;
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
                if(timeLupe%(TIMESKIP)==0)
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
    if(timeLupe%(TIMESKIP)==0)
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
