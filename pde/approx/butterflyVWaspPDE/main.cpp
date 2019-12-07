#include <string>

#include "numericaltrials.h"

#define SKIP_PRINT_UPDATE 40000
#define SKIP_FILE_SAVE 150
#define NUMBER_TIME_LOOP 30000000
#define MAX_NEWTON_STEPS 50
#define LEGENDRE_POLY_DEGREE 30
#define MAX_DELTA_NORM 0.0001
#define NUMBER_THREADS 15

int main()
{
    // Set up the temporal variables.
    double dt = 0.0001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    //double mu = 0.03; //0.01;
    double c  = 0.7;
    double g  = 2.0;
    double d  = 1.0;
    //double m  = 0.2;

    std::cout << "Starting" << std::endl;


//#define APPROXIMATE_MULTIPLE_M
#ifdef APPROXIMATE_MULTIPLE_M
    NumericalTrials::multipleApproximationsByM(
            mu,c,g,d,
            0.2,1.3,0.2,
            dt,NUMBER_TIME_LOOP,
            LEGENDRE_POLY_DEGREE,
            MAX_DELTA_NORM,MAX_NEWTON_STEPS,
            SKIP_PRINT_UPDATE,SKIP_FILE_SAVE,
            NUMBER_THREADS);
#else
    NumericalTrials trial;
    trial.approximateSystemTrackRepeating(
                0.02,0.04,9,
                c,g,d,
                0.05,1.2,110,
                dt,NUMBER_TIME_LOOP,
                LEGENDRE_POLY_DEGREE,
                MAX_DELTA_NORM,MAX_NEWTON_STEPS,
                -SKIP_PRINT_UPDATE,NUMBER_THREADS);
#endif

    std::cout << "Done" << std::endl;

    //delete trials;
    return(0);
}



