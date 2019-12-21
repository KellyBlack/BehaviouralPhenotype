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
    double dt = 0.00001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    double mu = 0.5;
    double c  = 2.8;
    double g  = 0.6;
    double d  = 0.1;
    //double m  = 0.2;

    std::cout << "Starting" << std::endl;


#define APPROXIMATE_MULTIPLE_M
#ifdef APPROXIMATE_MULTIPLE_M
    NumericalTrials::multipleApproximationsByM(
            mu,c,g,d,
            1.0,15.0,1.0,
            dt,NUMBER_TIME_LOOP,
            LEGENDRE_POLY_DEGREE,
            MAX_DELTA_NORM,MAX_NEWTON_STEPS,
            SKIP_PRINT_UPDATE,SKIP_FILE_SAVE,
            NUMBER_THREADS);
#else
    NumericalTrials trial;

    /*
    trial.approximateSystemHysteresis(
                mu,c,g,d,
                0.01,15.5,120,
                dt,NUMBER_TIME_LOOP,
                LEGENDRE_POLY_DEGREE,
                MAX_DELTA_NORM,MAX_NEWTON_STEPS,
                -SKIP_PRINT_UPDATE,
                false);

    */
    trial.approximateSystemTrackRepeating(
                0.01,0.15,4,
                c,g,d,
                0.01,15.5,120,
                dt,NUMBER_TIME_LOOP,
                LEGENDRE_POLY_DEGREE,
                MAX_DELTA_NORM,MAX_NEWTON_STEPS,
                -SKIP_PRINT_UPDATE,NUMBER_THREADS,
                false);

#endif

    std::cout << "Done" << std::endl;

    //delete trials;
    return(0);
}



