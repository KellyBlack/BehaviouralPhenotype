#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

#include "numericaltrials.h"

//#define OUTPUTFILE "approximation.csv"
#define BINARYOUTPUTFILE "approximation"
#define SKIP_PRINT_UPDATE 10000
#define SKIP_FILE_SAVE 150
#define NUMBER_TIME_LOOP 3000000
#define MAX_NEWTON_STEPS 50
#define LEGENDRE_POLY_DEGREE 30
#define MAX_DELTA_NORM 0.0001

int saveSeparateRunsByM(double mu,double c,double g,double d,
                        double lowM,double highM,double stepM,double dt,
                        NumericalTrials &trials);

int main()
{
    // Set up the temporal variables.
    double dt = 0.0001;

    NumericalTrials trials;

    // Define the default values of the parameters.
    double mu = 0.01;
    double c  = 0.7;
    double g  = 2.0;
    double d  = 1.0;
    //double m  = 0.2;

    std::cout << "Starting" << std::endl;
    saveSeparateRunsByM(mu,c,g,d,0.2,1.3,0.2,dt,trials);
    std::cout << "Done" << std::endl;

    return(0);
}


int saveSeparateRunsByM(double mu,double c,double g,double d,
                        double lowM,double highM,double stepM,double dt,
                        NumericalTrials &trials)
{

    for(double m=lowM;m<=highM;m+=stepM)
    {
        std::ostringstream filename("");
        filename << BINARYOUTPUTFILE
                 << "-m-" <<  std::setw(6) << std::fixed << std::setprecision(4) << std::setfill('0') << m
                 << "-mu-" << mu
                 << ".bin";
        std::cout << "Writing to " << filename.str() << std::endl;

        // Variables used to save the results of calculations into a
        // data file.
        std::fstream binFile (filename.str(), std::ios::out | std::ios::binary);

        if(trials.approximateSystem(
                    mu,c,g,d,m,
                    dt,NUMBER_TIME_LOOP,
                    LEGENDRE_POLY_DEGREE,
                    MAX_DELTA_NORM,MAX_NEWTON_STEPS,
                    binFile,
                    SKIP_PRINT_UPDATE,
                    SKIP_FILE_SAVE
          ) == 0)
            return(0);

        // Clean up the data file and close it
        binFile.close();

    }

    return(1);
}
