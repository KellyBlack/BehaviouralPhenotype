#include <fstream>
#include <sstream>
#include <math.h>

#include "rungakutta45.h"


// Constructor for the class. Sets the initial condition and time step.
RungaKutta45::RungaKutta45(double initialTimeStep, double *initialCondition)
{
    setDT(initialTimeStep);
    setCurrentTime(0.0);
    if(initialCondition!=nullptr)
    {
        state[0] = initialCondition[0];
        state[1] = initialCondition[1];
    }
    else
    {
        state[0] = 0.0;
        state[1] = 0.0;
    }

}

// Method to make a full approximation given a set of parameters.
// Also includes the initial condition and time span.
int RungaKutta45::approximation(double cValue,double gValue,double dValue,double mValue,double thetaValue,
                                double startTime,double endTime,double initialDt,
                                double *initialCond,double tolerance,
                                std::string filename,bool appendFile)
{
    // First set the values of the parameters.
    setC(cValue);
    setG(gValue);
    setD(dValue);
    setM(mValue);
    setTheta(thetaValue);
    setDT(initialDt);
    setCurrentTime(startTime);

    // Set the initial condition.
    state[0] = initialCond[0];
    state[1] = initialCond[1];

    // Open a file to write the results to.
    // Write the header for the file as well if this is a new file..
    std::fstream csvFile;
    if(appendFile)
    {
        csvFile.open(filename, std::ios::out | std::ios::app);
    }
    else
    {
        csvFile.open(filename, std::ios::out);
        csvFile << "which,mu,c,g,d,m,theta,time,maxWasp,minWasp,minButterfly,maxButterfly" << std::endl;
    }


    // Step through for the approximation. Stop when we hit the final time.q
    while(currentTime<endTime)
    {
        singleStep(tolerance);
    }

    // Life is good. End it now.
    csvFile.close();
    return(1);
}


// Method to calculate the value of the derivative of the state.
// Given a state and time calculate the rate of change.
void RungaKutta45::rateFunction(double,double *x,double *rate)
{
    rate[0] = dt*x[0];
    rate[1] = -dt*x[1];
}

// Method to make a single step of the Runga-Kutta Fehlberg scheme.
double RungaKutta45::singleStep(double tolerance)
{
    // Set up the vectors for the rate of change of the state
    // for the different stages.
    double rate0[2];
    double rate1[2];
    double rate2[2];
    double rate3[2];
    double rate4[2];
    double rate5[2];

    // Set up the vectors for the state at different stages.
    //double state0[2];
    double state1[2];
    double state2[2];
    double state3[2];
    double state4[2];
    double state5[2];

    // Set up the vectors that represent the increment for the
    // two approximations at the current times step.
    double updateFourth[2];
    double updateFifth[2];

    // Calculate the rate of change for the initial stage.
    //state0[0] = state[0];
    //state0[1] = state[1];
    rateFunction(currentTime,state,rate0);

    // Calculate the state and rate of change for the second stage.
    state1[0] = state[0] + 0.25*rate0[0];
    state1[1] = state[1] + 0.25*rate0[1];
    rateFunction(currentTime+0.25*dt,state1,rate1);

    // Calculate the state and rate of change for the third stage.
    state2[0] = state[0] + (3.0*rate0[0] + 9.0*rate1[0])/32.0;
    state2[1] = state[1] + (3.0*rate0[1] + 9.0*rate1[1])/32.0;
    rateFunction(currentTime+3.0/8.0*dt,state2,rate2);

    // Calculate the state and rate of change for the fourth stage.
    state3[0] = state[0] + (1932.0*rate0[0] + 7200.0*rate1[0] + 7296.0*rate2[0])/2197.0;
    state3[1] = state[1] + (1932.0*rate0[1] + 7200.0*rate1[1] + 7296.0*rate2[1])/2197.0;
    rateFunction(currentTime+12.0/13.0*dt,state3,rate3);

    // Calculate the state and rate of change for the fifth stage.
    state4[0] = state[0] + 439.0/216.0*rate0[0] - 8.0*rate1[0] + 3680.0/513.0*rate2[0] - 845.0/4104.0*rate3[0];
    state4[1] = state[1] + 439.0/216.0*rate0[1] - 8.0*rate1[1] + 3680.0/513.0*rate2[1] - 845.0/4104.0*rate3[1];
    rateFunction(currentTime+dt,state4,rate4);

    // Calculate the state and rate of change for the sixth stage.
    state5[0] = state[0] - 8.0/27.0*rate0[0] + 2.0*rate1[0] - 3544.0/2565.0*rate2[0] + 1859.0/4104.0*rate3[0] - 11.0/40.0*rate4[0];
    state5[1] = state[1] - 8.0/27.0*rate0[1] + 2.0*rate1[1] - 3544.0/2565.0*rate2[1] + 1859.0/4104.0*rate3[1] - 11.0/40.0*rate4[1];
    rateFunction(currentTime+0.5*dt,state5,rate5);

    // Calculate the update used to approximate the solution at the next time step.
    updateFourth[0] = 25.0/216.0*rate0[0] + 1408.0/2565.0*rate2[0] + 2197.0/4101.0*rate3[0] - 0.2*rate4[0];
    updateFourth[1] = 25.0/216.0*rate0[1] + 1408.0/2565.0*rate2[1] + 2197.0/4101.0*rate3[1] - 0.2*rate4[1];

    updateFifth[0] = 16.0/135.0*rate0[0] + 6656.0/12825.0*rate2[0] + 28561.0/56430.0*rate3[0] - 9.0/50.0*rate4[0] + 2.0/55.0*rate5[0];
    updateFifth[1] = 16.0/135.0*rate0[1] + 6656.0/12825.0*rate2[1] + 28561.0/56430.0*rate3[1] - 9.0/50.0*rate4[1] + 2.0/55.0*rate5[1];

    // Update the approximation for the next time step.
    state[0] += updateFourth[0];
    state[1] += updateFourth[1];
    currentTime += dt;

    // Now look at the difference between the two approximations
    // and update the value used for the time step.
    double s = pow(tolerance*dt*0.5/magnitudeDifference(updateFifth,updateFourth),0.25);
    dt = s*dt;

    return(dt);
}


// Method to approximate the magnitude of the difference between two
// state vectors. Returns a scalar for the magnitude of the difference.
double RungaKutta45::magnitudeDifference(double *one,double *two)
{
    double d1 = one[0]-two[0];
    double d2 = one[1]-two[1];
    return(sqrtf64(d1*d1+d2*d2));
}
