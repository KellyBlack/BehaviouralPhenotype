#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>

#include <sys/ipc.h>
#include <sys/msg.h>


#include "rungakutta45.h"


// Constructor for the class. Sets the initial condition and time step.
RungaKutta45::RungaKutta45(double initialTimeStep, double *initialCondition) :
    ApproximationBase::ApproximationBase()
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

long RungaKutta45::approximationByM(double cValue, double gValue, double dValue, double thetaValue,
                  double lowM, double highM, long numberM,
                  double startTime, double endTime, double initialDt, double minimumDT,
                  double *initialCond, double tolerance,
                  std::string filename, bool appendFile, int numberThreads)
{

    // Create an ID and then create a message queue that will be associated with the ID
    key_t msg_key = ftok("butterfly",1995);
    int msgID = msgget(msg_key,0666|IPC_CREAT);

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
        csvFile << "which,c,g,d,m,theta,time,maxWasp,minWasp,minButterfly,maxButterfly" << std::endl;
        //csvFile << "time,butterfly,wasp" << std::endl;
    }

    RungaKutta45::MaxMinBuffer msgValue;
    RungaKutta45::MessageInformation* processInformation;
    std::vector<RungaKutta45::MessageInformation*> processes;

    // Need to make sure there are no messages in the buffer left
    // over from previous runs that were prematurely terminated.
    // (I will not go into the details about the agony of learning this lesson.)
    while(msgrcv(msgID,&msgValue,sizeof(msgValue),0x4,IPC_NOWAIT)>0)
    {
        std::cout << "Previous message in the queue: " << msgValue.which << std::endl;
    }

    long lupe;
    long currentNumberProcesses = 0;
    for(lupe=0;lupe<numberM;++lupe)
    {
        double m = lowM + static_cast<double>(lupe)*(highM-lowM)/static_cast<double>(numberM-1);

        // Set up the data structure that will keep track of the values of the coefficient for this numerical trial.
        // Then start a new thread.
        RungaKutta45 *newApproximator = new RungaKutta45();
        RungaKutta45::MessageInformation *newProcess = new RungaKutta45::MessageInformation;
        newProcess->which = lupe;
        newProcess->theta = thetaValue;
        newProcess->c     = cValue;
        newProcess->g     = gValue;
        newProcess->d     = dValue;
        newProcess->m     = m;
        newProcess->trial = static_cast<ApproximationBase*>(newApproximator);
        newProcess->process =
                new std::thread(&RungaKutta45::approximation,newApproximator,
                                lupe,cValue,gValue,dValue,m,thetaValue,
                                startTime,endTime,initialDt,minimumDT,
                                initialCond,tolerance,msgID);
        processes.push_back(newProcess);
        std::cout << "starting " << lupe << "/" << m << ": " << newProcess->process << std::endl;


        if(++currentNumberProcesses>=numberThreads)
        {
            msgrcv(msgID,&msgValue,sizeof(msgValue),0x4,0);
            processInformation = findReturnedProcessParameters(msgValue,processes);

            csvFile << msgValue.which << ","
                    << msgValue.c << "," << msgValue.g << ","
                    << msgValue.d << ","  << msgValue.m << "," << msgValue.theta << "," << msgValue.endTime << ","
                    << msgValue.maxWasp << "," << msgValue.minWasp << ","
                    << msgValue.minButterfly << "," << msgValue.maxButterfly << std::endl;

            if(processInformation!= nullptr)
            {
                std::cout << msgValue.which << ","
                        << msgValue.c << "," << msgValue.g << ","
                        << msgValue.d << ","  << msgValue.m << "," << msgValue.theta << "," << msgValue.endTime << ","
                        << msgValue.maxWasp << "," << msgValue.minWasp << ","
                        << msgValue.minButterfly << "," << msgValue.maxButterfly << std::endl;
                delete processInformation->trial;
                processInformation->process->join();
                currentNumberProcesses -= 1;
            }
        }

    }

    std::cout << std::endl << std::endl << "All processes finished. " << currentNumberProcesses  << std::endl;

    // All of the process have been started. Now wait for them all to end and record the results.
    while(currentNumberProcesses>0)
    {
        std::cout << "Need to wait: " << currentNumberProcesses << std::endl;
        msgrcv(msgID,&msgValue,sizeof(msgValue),0x4,0);
        processInformation = findReturnedProcessParameters(msgValue,processes);

        csvFile << msgValue.which << ","
                << msgValue.c << "," << msgValue.g << ","
                << msgValue.d << ","  << msgValue.m << "," << msgValue.theta << "," << msgValue.endTime << ","
                << msgValue.maxWasp << "," << msgValue.minWasp << ","
                << msgValue.minButterfly << "," << msgValue.maxButterfly << std::endl;

        if(processInformation!= nullptr)
        {
            std::cout << msgValue.which << ","
                    << msgValue.c << "," << msgValue.g << ","
                    << msgValue.d << ","  << msgValue.m << "," << msgValue.theta << "," << msgValue.endTime << ","
                    << msgValue.maxWasp << "," << msgValue.minWasp << ","
                    << msgValue.minButterfly << "," << msgValue.maxButterfly << std::endl;
            processInformation->process->join();
            delete processInformation->trial;
            currentNumberProcesses -= 1;
        }
    }

    // Life is good. End it now.
    msgctl(msgID, IPC_RMID, nullptr);
    csvFile.close();


    return(lupe);

}

// Method to make a full approximation given a set of parameters.
// Also includes the initial condition and time span.
long RungaKutta45::approximation(long which,
                                 double cValue, double gValue, double dValue, double mValue, double thetaValue,
                                 double startTime, double endTime, double initialDt, double minimumDT,
                                 double *initialCond,double tolerance,
                                 int msgID)
{
    const long MAX_CHECK_CONSTANT = 4000;
    // First set the values of the parameters.
    setC(cValue);
    setG(gValue);
    setD(dValue);
    setM(mValue);
    setTheta(thetaValue);
    setDT(initialDt);
    setCurrentTime(startTime);
    setMinimumDT(minimumDT);

    // Set the initial condition.
    if(initialCond[0]<0.0)
    {
        state[0] = fabs(cValue*dValue/((gValue-dValue)*(1.0+mValue)))*1.0001;
        state[1] = fabs((1.0-state[0])*(cValue+state[0]*(1.0+mValue)));
    }
    else
    {
        state[0] = initialCond[0];
        state[1] = initialCond[1];
    }


    // Set up the variables used to determine the max and min of a cycle.
    double maxButterfliesDensity = state[0];
    double prevButterflyDensity = maxButterfliesDensity;
    double minButterfliesDensity = maxButterfliesDensity;
    int countButterflyIncreasing = 0;

    double prevWaspDensity = state[1];
    double maxWaspDensity = prevWaspDensity;
    double minWaspDensity = maxWaspDensity;
    int countWaspIncreasing = 0;
    int prevCycleClose = 0;

    double prevButterflyCheck = 0.0;
    double prevWaspCheck = 0.0;
    int prevValueClose = 0;

    // Step through for the approximation. Stop when we hit the final time.
    //unsigned long numberSteps = 0;
    while(currentTime<endTime)
    {
        singleStep(tolerance);
        /*
        if((numberSteps++)%2000==0)
        {
            std::cout << currentTime << "," << state[0] << "," << state[1] << std::endl;
            //csvFile << currentTime << "," << state[0] << "," << state[1] << std::endl;
        }
        */

        double currentButterflyDensity = state[0];
        double currentWaspDensity      = state[1];

        if((fabs(currentWaspDensity-prevWaspDensity)<1.0E-8)&&(fabs(currentButterflyDensity-prevButterflyDensity)<1.0E-8))
        {
            if(((fabs(prevWaspCheck-currentWaspDensity)>1.0E-8)||(fabs(prevButterflyCheck-currentButterflyDensity)>1.0E-8)))
            {
                prevValueClose = 0;
                prevWaspCheck = currentWaspDensity;
                prevButterflyCheck = currentButterflyDensity;
            }
            else
            {
                prevValueClose += 1;
            }

        }
        else
        {
            prevValueClose = 0;
            //prevWaspCheck = currentWaspDensity;
            //prevButterflyCheck = currentButterflyDensity;
        }

        if(currentButterflyDensity<prevButterflyDensity)
        {
            // The butterfly population is decreasing. If countIncreasing is
            // pos. that means that the butterfly count had a long trend of
            // increasing.
            if(countButterflyIncreasing-->0)
            {
                //std::cout << which << ": max butterfly " << prevButterflyDensity << "-" << maxButterfliesDensity << std::endl;
                if(fabs(maxButterfliesDensity-prevButterflyDensity)<1E-4)
                    prevCycleClose = (prevCycleClose | 0x1);
                else
                    prevCycleClose = 0;

                maxButterfliesDensity = prevButterflyDensity;
                countButterflyIncreasing = 0;
            }
        }
        else
        {
            // The butterfly population is increasing. If countIncreasing is
            // neg. that means that the butterfly count had a long trend of
            // decreasing.
            if(countButterflyIncreasing++ < 0)
            {
                //std::cout << which << ": min butterfly " << prevButterflyDensity << "-" << minButterfliesDensity << std::endl;
                if(fabs(minButterfliesDensity-prevButterflyDensity)<1E-4)
                    prevCycleClose = (prevCycleClose | 0x2);
                else
                    prevCycleClose = 0;

                minButterfliesDensity = prevButterflyDensity;
                countButterflyIncreasing = 0;
            }



        }
        prevButterflyDensity = currentButterflyDensity;


        if(currentWaspDensity<prevWaspDensity)
        {
            // The wasp population is decreasing. If countIncreasing is
            // pos. that means that the wasp count had a long trend of
            // increasing.
            if(countWaspIncreasing-->0)
            {
                //std::cout << which << ": max wasp " << prevWaspDensity << "-" << maxWaspDensity << std::endl;
                if(fabs(maxWaspDensity-prevWaspDensity)<1E-4)
                    prevCycleClose = (prevCycleClose | 0x4);
                else
                    prevCycleClose = 0;

                maxWaspDensity = prevWaspDensity;
                countWaspIncreasing = 0;
            }
        }
        else
        {
            // The wasp population is increasing. If countIncreasing is
            // neg. that means that the wasp count had a long trend of
            // decreasing.
            if(countWaspIncreasing++ < 0)
            {
                //std::cout << which << ": min wasp " << prevWaspDensity << "-" << minWaspDensity << std::endl;
                if(fabs(minWaspDensity-prevWaspDensity)<1E-4)
                    prevCycleClose = (prevCycleClose | 0x8);
                else
                    prevCycleClose = 0;

                minWaspDensity = prevWaspDensity;
                countWaspIncreasing = 0;
            }
        }
        prevWaspDensity = currentWaspDensity;

        if((prevCycleClose>=0xf)||(maxButterfliesDensity>10.0))
            break;
        else if (prevValueClose>MAX_CHECK_CONSTANT)
        {
            maxButterfliesDensity = currentButterflyDensity;
            minButterfliesDensity = currentButterflyDensity;
            maxWaspDensity = currentWaspDensity;
            minWaspDensity = currentWaspDensity;
            break;
        }

    }

    RungaKutta45::MaxMinBuffer values;
    values.which = which;
    values.mtype = 0x4;
    values.theta = theta;
    values.c     = c;
    values.g     = g;
    values.d     = d;
    values.m     = m;
    values.endTime = currentTime;
    values.maxWasp = maxWaspDensity;
    values.minWasp = minWaspDensity;
    values.maxButterfly = maxButterfliesDensity;
    values.minButterfly = minButterfliesDensity;
    msgsnd(msgID,&values,sizeof(values),0);


    return(1);
}


// Method to calculate the value of the derivative of the state.
// Given a state and time calculate the rate of change.
void RungaKutta45::rateFunction(double,double *x,double *rate)
{
    double p = 1.0+m*theta;
    rate[0] =  dt*x[0];
    rate[1] = -dt*x[1];

    rate[0] = dt*(
            p*x[0]*(1.0-x[0]) -
             x[1]*p*x[0]/(c+x[0]*p)
            );
    rate[1] = dt*(
         -d*x[1]+g*x[1]*p*x[0]/(c+x[0]*p)
            );
}

// Method to make a single step of the Runga-Kutta Fehlberg scheme.
double RungaKutta45::singleStep(double tolerance)
{
    // Set up the vectors for the rate of change of the state
    // for the different stages.
    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];
    double k5[2];
    double k6[2];

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
    rateFunction(currentTime,state,k1);

    // Calculate the state and rate of change for the second stage.
    state1[0] = state[0] + 0.25*k1[0];
    state1[1] = state[1] + 0.25*k1[1];
    rateFunction(currentTime+0.25*dt,state1,k2);

    // Calculate the state and rate of change for the third stage.
    state2[0] = state[0] + (3.0*k1[0] + 9.0*k2[0])/32.0;
    state2[1] = state[1] + (3.0*k1[1] + 9.0*k2[1])/32.0;
    rateFunction(currentTime+3.0/8.0*dt,state2,k3);

    // Calculate the state and rate of change for the fourth stage.
    state3[0] = state[0] + (1932.0*k1[0] - 7200.0*k2[0] + 7296.0*k3[0])/2197.0;
    state3[1] = state[1] + (1932.0*k1[1] - 7200.0*k2[1] + 7296.0*k3[1])/2197.0;
    rateFunction(currentTime+12.0/13.0*dt,state3,k4);

    // Calculate the state and rate of change for the fifth stage.
    state4[0] = state[0] + 439.0/216.0*k1[0] - 8.0*k2[0] + 3680.0/513.0*k3[0] - 845.0/4104.0*k4[0];
    state4[1] = state[1] + 439.0/216.0*k1[1] - 8.0*k2[1] + 3680.0/513.0*k3[1] - 845.0/4104.0*k4[1];
    rateFunction(currentTime+dt,state4,k5);

    // Calculate the state and rate of change for the sixth stage.
    state5[0] = state[0] - 8.0/27.0*k1[0] + 2.0*k2[0] - 3544.0/2565.0*k3[0] + 1859.0/4104.0*k4[0] - 11.0/40.0*k5[0];
    state5[1] = state[1] - 8.0/27.0*k1[1] + 2.0*k2[1] - 3544.0/2565.0*k3[1] + 1859.0/4104.0*k4[1] - 11.0/40.0*k5[1];
    rateFunction(currentTime+0.5*dt,state5,k6);

    // Calculate the update used to approximate the solution at the next time step.
    updateFourth[0] = 25.0/216.0*k1[0] + 1408.0/2565.0*k3[0] + 2197.0/4101.0*k4[0] - 0.2*k5[0];
    updateFourth[1] = 25.0/216.0*k1[1] + 1408.0/2565.0*k3[1] + 2197.0/4101.0*k4[1] - 0.2*k5[1];

    updateFifth[0] = 16.0/135.0*k1[0] + 6656.0/12825.0*k3[0] + 28561.0/56430.0*k4[0] - 9.0/50.0*k5[0] + 2.0/55.0*k6[0];
    updateFifth[1] = 16.0/135.0*k1[1] + 6656.0/12825.0*k3[1] + 28561.0/56430.0*k4[1] - 9.0/50.0*k5[1] + 2.0/55.0*k6[1];

    // Update the approximation for the next time step.
    state[0] += updateFourth[0];
    state[1] += updateFourth[1];
    currentTime += dt;

    // Now look at the difference between the two approximations
    // and update the value used for the time step.
    double s = pow(tolerance*0.5/magnitudeDifference(updateFifth,updateFourth),0.25);
    dt = s*dt;
    if(dt>minimumDT) dt = minimumDT;

    return(dt);
}


// Method to approximate the magnitude of the difference between two
// state vectors. Returns a scalar for the magnitude of the difference.
double RungaKutta45::magnitudeDifference(double *one,double *two)
{
    double d1 = one[0]-two[0];
    double d2 = one[1]-two[1];
    return(sqrt(d1*d1+d2*d2));
}
