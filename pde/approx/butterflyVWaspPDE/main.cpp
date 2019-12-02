#include <string>

#include "numericaltrials.h"

#define SKIP_PRINT_UPDATE 10000
#define SKIP_FILE_SAVE 150
#define NUMBER_TIME_LOOP 3000000
#define MAX_NEWTON_STEPS 50
#define LEGENDRE_POLY_DEGREE 30
#define MAX_DELTA_NORM 0.0001
#define NUMBER_THREADS 10

int main()
{
    // Set up the temporal variables.
    double dt = 0.0001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    double mu = 0.01;
    double c  = 0.7;
    double g  = 2.0;
    double d  = 1.0;
    //double m  = 0.2;

    std::cout << "Starting" << std::endl;

    NumericalTrials::multipleApproximationsByM(
            mu,c,g,d,
            0.2,1.3,0.2,
            dt,NUMBER_TIME_LOOP,
            LEGENDRE_POLY_DEGREE,
            MAX_DELTA_NORM,MAX_NEWTON_STEPS,
            SKIP_PRINT_UPDATE,SKIP_FILE_SAVE,
            NUMBER_THREADS);

    std::cout << "Done" << std::endl;

    //delete trials;
    return(0);
}



