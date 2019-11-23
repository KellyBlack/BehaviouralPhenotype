#include <iostream>
#include <fstream>
#include <cmath>

#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"
#include "butterflies.h"

//#define OUTPUTFILE "approximation.csv"
#define BINARYOUTPUTFILE "approximation.bin"
#define SKIP_PRINT_UPDATE 10000
#define SKIP_FILE_SAVE 75
#define NUMBER_TIME_LOOP 1000000
#define MAX_NEWTON_STEPS 50
#define LEGENDRE_POLY_DEGREE 30
#define MAX_DELTA_NORM 0.0001

int approximateSystem(double mu,double c,double g,double d,double a,
                      double dt,int maxTimeLupe,
                      Butterflies &theButterflies,
                      std::fstream &binFile)
{
    double t        = 0.0;
    int    timeLupe = 0;

    std::cout << "Pre-processing" << std::endl;
    theButterflies.initializeLegendreParams();
    theButterflies.setMu(mu);
    theButterflies.setC(c);
    theButterflies.setG(g);
    theButterflies.setD(d);
    theButterflies.setA(a);
    theButterflies.setDT(dt);

    theButterflies.initializeButterflies();
    theButterflies.writeParameters(binFile);
    theButterflies.writeBinaryHeader(binFile);

    // Start the time loop, and calculation an approximation at
    // each time step.
    double stepDeltaNorm = 0.0;
    for(timeLupe=0;(timeLupe<maxTimeLupe)&&(stepDeltaNorm<MAX_DELTA_NORM);++timeLupe)
    {
        t = static_cast<double>(timeLupe)*dt;

        if(timeLupe%(SKIP_PRINT_UPDATE)==0)
        {
            std::cout << "Calculating an approximation: " << timeLupe << " (" << t << ") ";
        }

        if(theButterflies.singleTimeStep(
                    MAX_DELTA_NORM,MAX_NEWTON_STEPS,timeLupe%(SKIP_PRINT_UPDATE)==0)
                < 0)
        {
            std::cout << std::endl << "Error - Newton's Method did not converge." << std::endl;
            return(0);
        }

        if(timeLupe%SKIP_FILE_SAVE==0)
        {
            theButterflies.writeBinaryCurrentApprox(t,binFile);
        }

    }


    return(1);
}

int main()
{
    // Set up the temporal variables.
    double dt       = 0.000185;
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
    std::fstream binFile (BINARYOUTPUTFILE, std::ios::out | std::ios::binary);

    double mu = 0.5;
    double c  = 0.7;
    double g  = 4.0;
    double d  = 2.0;
    double a  = 0.6;
    approximateSystem(mu,c,g,d,a,
                      dt,NUMBER_TIME_LOOP,
                      theButterflies,
                      binFile);

    // Clean up the data file and close it
#ifdef OUTPUTFILE
    resultsFile.close();
#endif
    binFile.close();

    std::cout << "Done" << std::endl;
    return(0);
}
