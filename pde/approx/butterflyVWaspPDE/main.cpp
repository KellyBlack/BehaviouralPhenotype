#include <string>

#include "numericaltrials.h"
#include "rungakutta45.h"

#define SKIP_PRINT_UPDATE 40000
#define SKIP_FILE_SAVE 150
#define NUMBER_TIME_LOOP 300000000
#define MAX_NEWTON_STEPS 50
#define LEGENDRE_POLY_DEGREE 30
#define MAX_DELTA_NORM 0.0001
#define NUMBER_THREADS 4

int main()
{
    // Set up the temporal variables.
    double dt = 0.00001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    double mu = 0.01;
    double c  = 1.1;
    double g  = 0.6;
    double d  = 0.1;
    double m  = 0.2;

    std::cout << "Starting" << std::endl;

//#define ODE_APPROXIMATION
#ifdef ODE_APPROXIMATION

    double theta = 1.0;
    double initialCond[2];
    RungaKutta45 odeApprox;

    initialCond[0] = c*d/((g-d)*(1.0+m))*0.95;
    initialCond[1] = (1.0-initialCond[0])*(c+initialCond[0]*(1.0+m));
    odeApprox.approximationByM(c,g,d,theta,
                               0.1,15.0,300,
                               0.0,500.0,dt,1.0E-5,
                               initialCond,1.0E-6,
                               "rk45_c1.1.csv",false,NUMBER_THREADS);

#endif

//#define APPROXIMATE_MULTIPLE_M
#ifdef APPROXIMATE_MULTIPLE_M
    NumericalTrials::multipleApproximationsByM(
            mu,c,g,d,
            //0.01,1.0,0.1,
            1.0,10.0,0.5,
            //10.5,12.0,0.5,
            //13.0,15.0,2.5,
            dt,NUMBER_TIME_LOOP,
            LEGENDRE_POLY_DEGREE,
            MAX_DELTA_NORM,MAX_NEWTON_STEPS,
            SKIP_PRINT_UPDATE,SKIP_FILE_SAVE,
            NUMBER_THREADS);
#endif


//#define APPROXIMATE_HYSTERESIS
#ifdef APPROXIMATE_HYSTERESIS
    NumericalTrials trial;

    trial.approximateSystemHysteresis(
                mu,c,g,d,
                //0.01,15.5,120,
                0.1,10.5,100,
                dt,NUMBER_TIME_LOOP,
                LEGENDRE_POLY_DEGREE,
                MAX_DELTA_NORM,MAX_NEWTON_STEPS,
                -SKIP_PRINT_UPDATE,
                true);

#endif

#define APPROXIMATE_MULTIPLE_M_MU
#ifdef APPROXIMATE_MULTIPLE_M_MU
    NumericalTrials trial;
    trial.approximateSystemTrackRepeating(
                mu,4.0*mu,4,
                c,g,d,
                //7.01,15.5,120,
                0.1,15.0,160,
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



