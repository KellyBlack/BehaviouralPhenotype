#ifndef MAINROUTINES_H
#define MAINROUTINES_H

#define SKIP_PRINT_UPDATE 100000000
#define SKIP_FILE_SAVE 1500
#define NUMBER_TIME_LOOP 3000000000
#define MAX_NEWTON_STEPS 50
#define LEGENDRE_POLY_DEGREE 30
#define MAX_DELTA_NORM 0.0001
#define NUMBER_THREADS 14

#include <string>
#include <sstream>

#include "numericaltrials.h"
#include "rungakutta45.h"

#include "numericaltrials.h"
//#include "util.h"
//#include "limInf.h"

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


void performManyApprpoximations_m()
{

    // Set up the temporal variables.
    double dt = 0.00001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    double mu = 0.01; // 0.0095
    double c  = 2.75;
    double g  = 0.6;
    double d  = 0.1;
    //double m  = 0.2;

    std::cout << "Starting" << std::endl;

    NumericalTrials::multipleApproximationsByM(
        mu,c,g,d,
        //0.01,1.0,0.1,
        6.5,13.1,3.25,
        //10.5,12.0,0.5,
        //13.0,15.0,2.5,
        dt,NUMBER_TIME_LOOP,
        LEGENDRE_POLY_DEGREE,
        MAX_DELTA_NORM,MAX_NEWTON_STEPS,
        SKIP_PRINT_UPDATE,SKIP_FILE_SAVE,
        NUMBER_THREADS);


}


void performManyApprpoximations_m_c()
{

    // Set up the temporal variables.
    double dt = 0.00001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    double mu = 0.01; // 0.0095
    double c  = 2.75;
    double g  = 0.6;
    double d  = 0.1;
    double m  = 0.2;

    std::ostringstream filename("");

    std::cout << "Starting" << std::endl;

    /*
    NumericalTrials::multipleApproximationsByMandC(
                mu,g,d,
                0.1,10.0,0.1,
                0.1,50.0,0.1,
                //2.0,2.1,0.1,
                //17.0,22.01,0.1,
                dt,NUMBER_TIME_LOOP,
                LEGENDRE_POLY_DEGREE,
                MAX_DELTA_NORM,MAX_NEWTON_STEPS,
                SKIP_PRINT_UPDATE,SKIP_FILE_SAVE,
                NUMBER_THREADS,
                "/tmp/results.csv");
    */

    struct remoteProcess {
        std::thread process;
        std::atomic<bool> running;
        NumericalTrials* trial;
    };

    NumericalTrials *trial    = new NumericalTrials();
    remoteProcess* newProcess = new remoteProcess;
    newProcess->running = true;
    newProcess->trial = trial;
    c = 0.7;
    m = 10.0;
    //trial->approximateSystem(
    //            mu,c,g,d,m,
    //            dt,NUMBER_TIME_LOOP,
    //            LEGENDRE_POLY_DEGREE,
    //            MAX_DELTA_NORM,MAX_NEWTON_STEPS,
    //            "/tmp/testing.bin",
    //            SKIP_PRINT_UPDATE,SKIP_FILE_SAVE);
    trial->approximateSystemCheckOscillation(
        mu,c,g,d,m,
        dt,NUMBER_TIME_LOOP,
        LEGENDRE_POLY_DEGREE,
        MAX_DELTA_NORM,MAX_NEWTON_STEPS,
        "/tmp/testing.bin",
        SKIP_PRINT_UPDATE,SKIP_FILE_SAVE,
        nullptr,&(newProcess->running)
        );
    delete trial;


}

void makeOneApproximation()
{
    // Set up the temporal variables.
    double dt = 0.00001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    double mu = 0.01; // 0.0095
    double c  = 2.75;
    double g  = 0.6;
    double d  = 0.1;
    double m  = 0.2;

    std::cout << "Starting" << std::endl;

    NumericalTrials trial;
    c = 2.0;  // 2.0; // 8.0;
    m = 10.0; //18.0;  // 5.0;
    d = 0.1;
    g = 0.6;
    mu = 0.01;
    trial.approximateSystem(
        mu,c,g,d, m,
        dt,NUMBER_TIME_LOOP,
        LEGENDRE_POLY_DEGREE,
        MAX_DELTA_NORM,MAX_NEWTON_STEPS,
        "/tmp/trialApprox.bin", //approximation-c-2.0-m-18.0-mu-0.01.bin",
        SKIP_PRINT_UPDATE,SKIP_FILE_SAVE);

    /*
    struct remoteProcess {
        std::thread process;
        std::atomic<bool> running;
        NumericalTrials* trial;
    };
    NumericalTrials *trial =  new NumericalTrials();
    remoteProcess* newProcess = new remoteProcess;
    newProcess->running = true;
    newProcess->trial = trial;
    trial->approximateSystemCheckOscillation(
            mu, c, g, d, m,
            dt, NUMBER_TIME_LOOP,
            LEGENDRE_POLY_DEGREE,
            MAX_DELTA_NORM,MAX_NEWTON_STEPS,
            "/tmp/trial.dat",
            SKIP_PRINT_UPDATE,SKIP_FILE_SAVE,
            nullptr,
            &(newProcess->running));
    */

}

void checkHysteresis_by_m()
{
    // Set up the temporal variables.
    double dt = 0.00001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    double mu = 0.01; // 0.0095
    double c  = 2.75;
    double g  = 0.6;
    double d  = 0.1;
    //double m  = 0.2;

    std::cout << "Starting" << std::endl;

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


}


void performManyApproximations_by_m_my()
{
    // Set up the temporal variables.
    double dt = 0.00001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    double mu = 0.01; // 0.0095
    double c  = 2.75;
    double g  = 0.6;
    double d  = 0.1;
    //double m  = 0.2;

    std::cout << "Starting" << std::endl;

    std::ostringstream filename("");
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

}

#endif // MAINROUTINES_H
