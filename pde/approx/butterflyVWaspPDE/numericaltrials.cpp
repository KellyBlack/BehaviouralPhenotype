#include <fstream>
#include <vector>
#include <thread>
#include <string>
#include <sstream>

#include <sys/ipc.h>
#include <sys/msg.h>


#include "numericaltrials.h"
#include "util.h"

//#define OUTPUTFILE "approximation.csv"
#define BINARYOUTPUTFILE "approximation"


// Constructor for the numerical trials class.
NumericalTrials::NumericalTrials()
{

}

void NumericalTrials::multipleApproximationsByM(
        double mu, double c, double g, double d,
        double lowM, double highM, double stepM,
        double dt, int maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint,int skipFileSave,
        int numberThreads)
{

    int num = 0;
    double m = lowM;
    std::vector<std::thread> processes;
    std::vector<NumericalTrials*> trials;
    while(m<=highM)
    {
        processes.clear();
        trials.clear();
        num = 0;
        while((m<=highM)&&(num<numberThreads))
        {
            std::ostringstream filename("");
            filename << BINARYOUTPUTFILE
                     << "-m-" <<  std::setw(6) << std::fixed << std::setprecision(4) << std::setfill('0') << m
                     << "-mu-" << mu
                     << ".bin";
            std::cout << "Writing to " << filename.str() << std::endl;

            NumericalTrials *trial = new NumericalTrials();
            trials.push_back(trial);
            processes.push_back(
                        std::thread(&NumericalTrials::approximateSystem,trial,
                                    mu,c,g,d,m,
                                    dt,maxTimeLupe,
                                    legendrePolyDegree,
                                    maxDeltaNorm,maxNewtonSteps,
                                    filename.str(),
                                    skipPrint,
                                    skipFileSave)
                        );

            std::cout << "Moving along " << m << std::endl;
            m += stepM;
            num += 1;
        }

        std::vector<std::thread>::iterator eachProcess;
        for(eachProcess=processes.begin();eachProcess!=processes.end();++eachProcess)
        {
            eachProcess->join();
        }

        std::vector<NumericalTrials*>::iterator eachTrial;
        for(eachTrial=trials.begin();eachTrial!=trials.end();++eachTrial)
        {
            delete *eachTrial;
        }
        //th.join();
    }


}

int NumericalTrials::approximateSystem(double mu, double c, double g, double d, double m,
        double dt, int maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        std::string filename,
        int skipPrint, int skipFileSave)
{
    double t        = 0.0;
    int    timeLupe = 0;

    // Variables used to save the results of calculations into a
    // data file.
    std::fstream binFile (filename, std::ios::out | std::ios::binary);

    std::cout << "Pre-processing" << std::endl;
    int N = legendrePolyDegree;
    Butterflies *theButterflies = new Butterflies(N,N+2);
    theButterflies->initializeLegendreParams();
    theButterflies->setMu(mu);
    theButterflies->setC(c);
    theButterflies->setG(g);
    theButterflies->setD(d);
    theButterflies->setM(m);
    theButterflies->setDT(dt);

    theButterflies->initializeButterfliesGaussian(0.0,mu);
    theButterflies->writeParameters(binFile);
    theButterflies->writeBinaryHeader(binFile);

    // Start the time loop, and calculation an approximation at
    // each time step.
    for(timeLupe=0;timeLupe<maxTimeLupe;++timeLupe)
    {
        t = static_cast<double>(timeLupe)*dt;

        if(timeLupe%(skipPrint)==0)
        {
            std::cout << "Calculating an approximation: "
                         << std::fixed
                         << std::setw(8)
                         << std::setprecision(4)
                         << timeLupe << " (" << t << ") ";
        }

        if(
           theButterflies->singleTimeStep(maxDeltaNorm,maxNewtonSteps,timeLupe%(skipPrint)==0)
                < 0)
        {
            std::cout << std::endl << "Error - Newton's Method did not converge." << std::endl;
            return(0);
        }

        if(timeLupe%skipFileSave==0)
        {
            theButterflies->writeBinaryCurrentApprox(t,binFile);
        }

    }

    // Clean up the data file and close it
    binFile.close();
    delete theButterflies;

    return(1);
}

int NumericalTrials::approximateSystemTrackRepeating(
        double mu, double c, double g, double d, double m,
        double dt, int maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint)
{

    // Create an ID and then create a message queue that will be associated with the ID
    key_t msg_key = ftok("butterfly",1995);
    int msgID = msgget(msg_key,0666|IPC_CREAT);

    // Set up the vector that will be used to keep track of the processes and
    // returned information associated with each process.
    struct MessageInformation
    {
        unsigned long which;
        double mu;
        double c;
        double g;
        double d;
        double m;
        double time;
        double maxButterfly;
        double minButterfly;
        double maxWasp;
        double minWasp;
        std::thread *process;
        NumericalTrials *trial;
    };

    std::vector<MessageInformation*> processes;
    unsigned long lupe;
    for(lupe=0;lupe<8;++lupe)
    {
        std::cout << "Starting process " << lupe << std::endl;
        m = 0.01 + static_cast<double>(lupe)*0.05;
        MessageInformation *newProcess = new MessageInformation;
        newProcess->which = lupe;
        newProcess->mu = mu;
        newProcess->c = c;
        newProcess->g = g;
        newProcess->d = d;
        newProcess->m = m;
        newProcess->trial = new NumericalTrials;
        newProcess->process = new std::thread(
                        &NumericalTrials::approximateSystemQuietResponse,newProcess->trial,
                        mu,c,g,d,m,dt,maxTimeLupe/10,
                        legendrePolyDegree,maxDeltaNorm,maxNewtonSteps,skipPrint,msgID,lupe
                        );
        processes.push_back(newProcess);
    }

    MaxMinBuffer msgValue;
    std::vector<MessageInformation*>::iterator eachProcess;
    std::cout << "waiting on " << processes.size() << " processes." << std::endl;
    for(lupe=0;lupe<processes.size();++lupe)
    {
        std::cout << "waiting on response " << lupe << std::endl;
        msgrcv(msgID,&msgValue,sizeof(msgValue),2,0);
        std::cout << "Value: " << msgValue.which << " " << msgValue.m << " "
                  << msgValue.maxWasp << " " << msgValue.minWasp << " "
                  << msgValue.minButterfly << " " << msgValue.maxButterfly << std::endl;

        // need to find this thread and join it.
        for(eachProcess=processes.begin();(eachProcess!=processes.end());++eachProcess)
        {
            if((*eachProcess)->which == msgValue.which)
            {
                (*eachProcess)->mu   = msgValue.mu;
                (*eachProcess)->c    = msgValue.c;
                (*eachProcess)->g    = msgValue.g;
                (*eachProcess)->d    = msgValue.d;
                (*eachProcess)->m    = msgValue.m;
                (*eachProcess)->time = msgValue.endTime;
                (*eachProcess)->maxButterfly = msgValue.maxButterfly;
                (*eachProcess)->minButterfly = msgValue.minButterfly;
                (*eachProcess)->maxWasp      = msgValue.maxWasp;
                (*eachProcess)->minWasp      = msgValue.minWasp;

                std::cout << "Found the process " << (*eachProcess)->which << std::endl;
                delete (*eachProcess)->trial;
                (*eachProcess)->process->join();
                break;
            }
        }
        std::cout << "heard response " << lupe << " (" << processes.size() << ")" << std::endl;
    }

    std::cout << std::endl << std::endl << "All processes finished." << std::endl;
    for(eachProcess=processes.begin();eachProcess!=processes.end();++eachProcess)
    {
        //eachProcess->join();
        //delete *eachTrial;
        std::cout << "Value: " << (*eachProcess)->which << " "
                  << (*eachProcess)->mu << "," << (*eachProcess)->c << "," << (*eachProcess)->g << ","
                  << (*eachProcess)->d << "," << (*eachProcess)->m << "," << (*eachProcess)->time << ","
                  << (*eachProcess)->maxWasp << "," << (*eachProcess)->minWasp << ","
                  << (*eachProcess)->minButterfly << "," << (*eachProcess)->maxButterfly << std::endl;
    }


    /*
    NumericalTrials t1;
    std::thread p1(
        &NumericalTrials::approximateSystemQuietResponse,&t1,mu,c,g,d,m,dt,maxTimeLupe/10,
        legendrePolyDegree,maxDeltaNorm,maxNewtonSteps,skipPrint,msgID,1
                );

    NumericalTrials t2;
    std::thread p2(
        &NumericalTrials::approximateSystemQuietResponse,&t2,mu,c,g,d,m+0.1,dt,maxTimeLupe/10,
        legendrePolyDegree,maxDeltaNorm,maxNewtonSteps,skipPrint,msgID,2
                );

    MaxMinBuffer msgValue;
    msgrcv(msgID,&msgValue,sizeof(msgValue),2,0);

    std::cout << "Value: " << msgValue.which << " " << msgValue.maxWasp << " " << msgValue.minWasp << " "
              << msgValue.minButterfly << " " << msgValue.maxButterfly << std::endl;

    msgrcv(msgID,&msgValue,sizeof(msgValue),2,0);

    std::cout << "Value: " << msgValue.which << " " << msgValue.maxWasp << " " << msgValue.minWasp << " "
              << msgValue.minButterfly << " " << msgValue.maxButterfly << std::endl;


    p1.join();
    p2.join();
    */

    msgctl(msgID, IPC_RMID, nullptr);
    return(1);
}

int NumericalTrials::approximateSystemQuietResponse(
        double mu, double c, double g, double d, double m,
        double dt, int maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint,int msgID,unsigned long which)
{
    ArrayUtils<double> arrays;
    double *maxButterflyProfile = arrays.onetensor(legendrePolyDegree+2);
    double t        = 0.0;
    int    timeLupe = 0;


    //std::cout << "Pre-processing" << std::endl;
    int N = legendrePolyDegree;
    Butterflies *theButterflies = new Butterflies(N,N+2);
    theButterflies->initializeLegendreParams();
    theButterflies->setMu(mu);
    theButterflies->setC(c);
    theButterflies->setG(g);
    theButterflies->setD(d);
    theButterflies->setM(m);
    theButterflies->setDT(dt);

    theButterflies->initializeButterflies();
    theButterflies->copyCurrentState(maxButterflyProfile);
    double maxButterfliesDensity = theButterflies->totalButterflyPopulation();
    double prevButterflyDensity = maxButterfliesDensity;
    double minButterfliesDensity = maxButterfliesDensity;
    int countButterflyIncreasing = 0;

    double prevWaspDensity = theButterflies->waspPopulation();
    double maxWaspDensity = prevWaspDensity;
    double minWaspDensity = maxWaspDensity;
    int countWaspIncreasing = 0;


    // Start the time loop, and calculation an approximation at
    // each time step.
    for(timeLupe=0;timeLupe<maxTimeLupe;++timeLupe)
    {
        t = static_cast<double>(timeLupe)*dt;

        if((skipPrint>0) && (timeLupe%(skipPrint)==0))
        {
            std::cout << "Calculating an approximation: "
                         << std::fixed
                         << std::setw(8)
                         << std::setprecision(4)
                         << timeLupe << " (" << t << "-" << which <<  ") ";
        }

        if(
           theButterflies->singleTimeStep(maxDeltaNorm,maxNewtonSteps,(skipPrint>0)&&(timeLupe%(skipPrint)==0))
                < 0)
        {
            std::cout << std::endl << "Error - Newton's Method did not converge." << std::endl;
            return(0);
        }

        double currentButterflyDensity = theButterflies->totalButterflyPopulation();
        if(currentButterflyDensity<prevButterflyDensity)
        {
            // The butterfly population is decreasing. If countIncreasing is
            // pos. that means that the butterfly count had a long trend of
            // increasing.
            if(countButterflyIncreasing-->0)
            {
                //std::cout << which << ": max butterfly " << prevButterflyDensity << "-" << maxButterfliesDensity << std::endl;
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
                minButterfliesDensity = prevButterflyDensity;
                countButterflyIncreasing = 0;
            }



        }
        prevButterflyDensity = currentButterflyDensity;


        double currentWaspDensity = theButterflies->waspPopulation();
        if(currentWaspDensity<prevWaspDensity)
        {
            // The wasp population is decreasing. If countIncreasing is
            // pos. that means that the wasp count had a long trend of
            // increasing.
            if(countWaspIncreasing-->0)
            {
                //std::cout << which << ": max wasp " << prevWaspDensity << "-" << maxWaspDensity << std::endl;
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
                minWaspDensity = prevWaspDensity;
                countWaspIncreasing = 0;
            }
        }
        prevWaspDensity = currentWaspDensity;

    }

    MaxMinBuffer values;
    values.mtype = 2;
    values.mu   = mu;
    values.c    = c;
    values.g    = g;
    values.d    = d;
    values.m    = m;
    values.endTime = t;
    values.maxWasp = maxWaspDensity;
    values.minWasp = minWaspDensity;
    values.maxButterfly = maxButterfliesDensity;
    values.minButterfly = minButterfliesDensity;
    values.which = which;
    msgsnd(msgID,&values,sizeof(values),0);
    //std::cout << "DONE " << which << std::endl;

    // Clean up the allocated space
    delete theButterflies;
    arrays.delonetensor(maxButterflyProfile);

    return(1);

}

