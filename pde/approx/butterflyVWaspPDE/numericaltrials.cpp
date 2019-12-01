#include <fstream>

#include "numericaltrials.h"

NumericalTrials::NumericalTrials()
{

}

int NumericalTrials::approximateSystem(double mu, double c, double g, double d, double m,
        double dt, int maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        std::string filename,
        int skipPrint, int skipFileSave)
{
    double t        = 0.0;
    int    timeLupe = 0;

    // Variables used to save the results of calculations into a
    // data file.
    std::fstream binFile (filename, std::ios::out | std::ios::binary);

    std::cout << "Pre-processing" << std::endl;
    int N = legendrePolyDegree;
    Butterflies *theButterflies = new Butterflies(N,N+2);
    theButterflies->initializeLegendreParams();
    theButterflies->setMu(mu);
    theButterflies->setC(c);
    theButterflies->setG(g);
    theButterflies->setD(d);
    theButterflies->setM(m);
    theButterflies->setDT(dt);

    theButterflies->initializeButterflies();
    theButterflies->writeParameters(binFile);
    theButterflies->writeBinaryHeader(binFile);

    // Start the time loop, and calculation an approximation at
    // each time step.
    for(timeLupe=0;timeLupe<maxTimeLupe;++timeLupe)
    {
        t = static_cast<double>(timeLupe)*dt;

        if(timeLupe%(skipPrint)==0)
        {
            std::cout << "Calculating an approximation: "
                         << std::fixed
                         << std::setw(8)
                         << std::setprecision(4)
                         << timeLupe << " (" << t << ") ";
        }

        if(
           theButterflies->singleTimeStep(maxDeltaNorm,maxNewtonSteps,timeLupe%(skipPrint)==0)
                < 0)
        {
            std::cout << std::endl << "Error - Newton's Method did not converge." << std::endl;
            return(0);
        }

        if(timeLupe%skipFileSave==0)
        {
            theButterflies->writeBinaryCurrentApprox(t,binFile);
        }

    }

    // Clean up the data file and close it
    binFile.close();
    delete theButterflies;

    return(1);
}
