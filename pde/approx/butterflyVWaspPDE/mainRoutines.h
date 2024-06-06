#ifndef MAINROUTINES_H
#define MAINROUTINES_H

#define NUMBER_THREADS 14

#include <string>
#include <sstream>

#include "numericaltrials.h"
#include "rungakutta45.h"

#include "numericaltrials.h"
#include "util.h"
#include "limInf.h"

void odeApproximation()
{
    // Set up the temporal variables.
    double dt = 0.00001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    double c  = 2.75;
    double g  = 0.6;
    double d  = 0.1;
    double m  = 0.2;

    std::ostringstream filename("");

    std::cout << "Starting" << std::endl;


    filename.str("");
    filename.clear();
    filename << "/tmp/rk45_c-"
             <<  std::setw(6) << std::fixed << std::setprecision(4) << std::setfill('0') << c
             << ".csv";
    std::cout << "Writing to " << filename.str() << std::endl;

    double theta = 1.0;
    double initialCond[2];
    RungaKutta45 odeApprox;

    initialCond[0] = c*d/((g-d)*(1.0+m))*0.95;
    initialCond[1] = (1.0-initialCond[0])*(c+initialCond[0]*(1.0+m));
    odeApprox.approximationByM(c,g,d,theta,
                               0.1,15.0,300,
                               0.0,500.0,dt,1.0E-5,
                               initialCond,1.0E-6,
                               filename.str(),false,NUMBER_THREADS);

}

#endif // MAINROUTINES_H
