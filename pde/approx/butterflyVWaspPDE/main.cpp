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

     std::cout << "Starting" << std::endl;

    odeApproximation();
    //performManyApprpoximations_m();
    //performManyApprpoximations_m_c();
    //makeOneApproximation();
    //checkHysteresis_by_m();
    performManyApproximations_by_m_mu();

    //determineSteadyState();

    std::cout << "Done" << std::endl;

    //delete trials;
    return(0);
}



