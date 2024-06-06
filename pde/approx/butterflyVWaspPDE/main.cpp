//#include <string>
//#include <sstream>

//#include "numericaltrials.h"
//#include "rungakutta45.h"

//#include "numericaltrials.h"
//#include "util.h"
//#include "limInf.h"

#include "mainRoutines.h"



int main()
{

    //odeApproximation();
    //performManyApprpoximations_m();
    //performManyApprpoximations_m_c();
    //makeOneApproximation();
    checkHysteresis_by_m();

//#define APPROXIMATE_MULTIPLE_M_MU
#ifdef APPROXIMATE_MULTIPLE_M_MU

    filename.str("");
    filename.clear();
    filename << "/tmp/changingMResults_c="
             <<  std::setw(6) << std::fixed << std::setprecision(4) << std::setfill('0') << c
             << ".csv";
    std::cout << "Writing to " << filename.str() << std::endl;

    NumericalTrials trial;
    trial.approximateSystemTrackRepeating(
                mu,2.5*mu,6,
                c,g,d,
                //7.01,15.5,120,
                0.1,15.0,160,
                dt,NUMBER_TIME_LOOP,
                LEGENDRE_POLY_DEGREE,
                MAX_DELTA_NORM,MAX_NEWTON_STEPS,
                -SKIP_PRINT_UPDATE,NUMBER_THREADS,
                false,
                filename.str());
#endif

    std::cout << "Done" << std::endl;

    //delete trials;
    return(0);
}



