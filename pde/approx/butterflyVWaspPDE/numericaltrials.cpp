#include <fstream>
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

std::vector<NumericalTrials::MessageInformation*>::iterator NumericalTrials::findReturnedProcessParameters(
        NumericalTrials::MaxMinBuffer msgValue,
        std::vector<NumericalTrials::MessageInformation*> processes)
{

    std::vector<NumericalTrials::MessageInformation*>::iterator eachProcess;

    // need to find this thread and join it.
    for(eachProcess=processes.begin();(eachProcess!=processes.end());++eachProcess)
    {
        if((*eachProcess)->which == msgValue.which)
        {
            // This is the thread from which this process sprang forth.
            // record the values that were passed and clean up the mess.
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
        }
    }

    return(eachProcess);
}


int NumericalTrials::approximateSystemTrackRepeating(
        double muLow, double muHigh, int numberMu,
        double c, double g, double d,
        double mLow, double mHigh, int numberM,
        double dt, int maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint, int numProcesses,
        bool appendFile)
{

    // Create an ID and then create a message queue that will be associated with the ID
    key_t msg_key = ftok("butterfly",1995);
    int msgID = msgget(msg_key,0666|IPC_CREAT);

    // Open a file to write the results to.  Write the header for the file as well if this is a new file..
    std::fstream csvFile;
    if(appendFile)
    {
        csvFile.open("changingMResults.csv", std::ios::out | std::ios::app);
    }
    else
    {
        csvFile.open("changingMResults.csv", std::ios::out);
        csvFile << "which,mu,c,g,d,m,time,maxWasp,minWasp,minButterfly,maxButterfly" << std::endl;
    }

    double currentM = mLow;
    double currentDiffusion = muLow;

    // Set up the vector used to keep track of all the information that the
    // remote processes create.
    MaxMinBuffer msgValue;
    std::vector<NumericalTrials::MessageInformation*>::iterator eachProcess;
    std::vector<NumericalTrials::MessageInformation*> processes;

    // Need to make sure there are no messages in the buffer left
    // over from previous runs that were prematurely terminated.
    // (I will not go into the details about the agony of learning this lesson.)
    while(msgrcv(msgID,&msgValue,sizeof(msgValue),0,IPC_NOWAIT)>0)
    {
        std::cout << "Previous message in the queue: " << msgValue.which << std::endl;
    }

    // Loop through all possible values of the diffusion and value of m.
    unsigned long mLupe;
    unsigned long diffusionLupe;
    unsigned long currentNumberProcesses = 0;
    unsigned long totalRuns = 0;
    for(mLupe=0;mLupe<static_cast<unsigned long>(numberM);++mLupe)
        for(diffusionLupe=0;diffusionLupe<static_cast<unsigned long>(numberMu);++diffusionLupe)
        {
            // Calculate the value of m and the diffusion coefficient.
            currentM = mLow + static_cast<double>(mLupe)*(mHigh-mLow)/static_cast<double>(numberM-1);
            currentDiffusion = muLow + static_cast<double>(diffusionLupe)*(muHigh-muLow)/static_cast<unsigned long>(numberMu-1);

            // Set up the data structure that will keep track of the values of the coefficient for this numerical trial.
            // Then start a new thread.
            NumericalTrials::MessageInformation *newProcess = new NumericalTrials::MessageInformation;
            newProcess->which = totalRuns;
            newProcess->mu    = currentDiffusion;
            newProcess->c     = c;
            newProcess->g     = g;
            newProcess->d     = d;
            newProcess->m     = currentM;
            newProcess->trial = new NumericalTrials;
            newProcess->process = new std::thread(
                            &NumericalTrials::approximateSystemQuietResponse,newProcess->trial,
                            currentDiffusion,c,g,d,currentM,dt,maxTimeLupe,
                            legendrePolyDegree,maxDeltaNorm,maxNewtonSteps,skipPrint,msgID,totalRuns
                            );
            processes.push_back(newProcess);
            std::cout << "starting " << currentDiffusion << "/" << currentM << ": " << totalRuns << std::endl;

            totalRuns += 1;
            if(++currentNumberProcesses>=static_cast<unsigned long>(numProcesses))
            {
                // There are too many processes running. Need to wait for one to stop
                // before starting a new process.
                msgrcv(msgID,&msgValue,sizeof(msgValue),0,0);

                // One just ended. Figure out the values that were sent and record the information from the run.
                std::cout << msgValue.which << "--"
                          << msgValue.mu << "," << msgValue.c << "," << msgValue.g << ","
                          << msgValue.d << "," << msgValue.m << "," << msgValue.endTime << ","
                          << msgValue.maxWasp << "," << msgValue.minWasp << ","
                          << msgValue.minButterfly << "," << msgValue.maxButterfly << std::endl;

                // need to find this thread and join it.
                for(eachProcess=processes.begin();(eachProcess!=processes.end());++eachProcess)
                {
                    if((*eachProcess)->which == msgValue.which)
                    {
                        // This is the thread from which this process sprang forth.
                        // record the values that were passed and clean up the mess.
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

                        csvFile << (*eachProcess)->which << ","
                                << (*eachProcess)->mu << "," << (*eachProcess)->c << "," << (*eachProcess)->g << ","
                                << (*eachProcess)->d << "," << (*eachProcess)->m << "," << (*eachProcess)->time << ","
                                << (*eachProcess)->maxWasp << "," << (*eachProcess)->minWasp << ","
                                << (*eachProcess)->minButterfly << "," << (*eachProcess)->maxButterfly << std::endl;

                        delete (*eachProcess)->trial;
                        (*eachProcess)->process->join();
                        currentNumberProcesses -= 1;
                        break;
                    }
                }

            }
        }
    std::cout << std::endl << std::endl << "All processes finished. " << currentNumberProcesses  << std::endl;

    // All of the process have been started. Now wait for them all to end and record the results.
    while(currentNumberProcesses>0)
    {
        std::cout << "Need to wait: " << currentNumberProcesses << std::endl;
        msgrcv(msgID,&msgValue,sizeof(msgValue),0,0);
        std::cout << msgValue.which << " "
                  << msgValue.mu << "," << msgValue.c << "," << msgValue.g << ","
                  << msgValue.d << "," << msgValue.m << "," << msgValue.endTime << ","
                  << msgValue.maxWasp << "," << msgValue.minWasp << ","
                  << msgValue.minButterfly << "," << msgValue.maxButterfly << std::endl;

        // need to find this thread and join it.
        for(eachProcess=processes.begin();(eachProcess!=processes.end());++eachProcess)
        {
            if((*eachProcess)->which == msgValue.which)
            {
                // This is the thread from which this process sprang forth.
                // record the values that were passed and clean up the mess.
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

                csvFile << (*eachProcess)->which << ","
                        << (*eachProcess)->mu << "," << (*eachProcess)->c << "," << (*eachProcess)->g << ","
                        << (*eachProcess)->d << "," << (*eachProcess)->m << "," << (*eachProcess)->time << ","
                        << (*eachProcess)->maxWasp << "," << (*eachProcess)->minWasp << ","
                        << (*eachProcess)->minButterfly << "," << (*eachProcess)->maxButterfly << std::endl;

                delete (*eachProcess)->trial;
                (*eachProcess)->process->join();
                currentNumberProcesses -= 1;
                break;
            }
        }
    }

    for(eachProcess=processes.begin();eachProcess!=processes.end();++eachProcess)
    {
        std::cout << (*eachProcess)->which << ","
                  << (*eachProcess)->mu << "," << (*eachProcess)->c << "," << (*eachProcess)->g << ","
                  << (*eachProcess)->d << "," << (*eachProcess)->m << "," << (*eachProcess)->time << ","
                  << (*eachProcess)->maxWasp << "," << (*eachProcess)->minWasp << ","
                  << (*eachProcess)->minButterfly << "," << (*eachProcess)->maxButterfly << std::endl;
    }

    // Clean up all of the information and close the communication channel.
    processes.clear();
    csvFile.close();
    msgctl(msgID, IPC_RMID, nullptr);

    // Our work here is done. Live long and kick butt.
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


    //std::cout << "Pre-processing: " << which << std::endl;
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
    //theButterflies->initializeButterfliesGaussian(1.0,mu*0.5);
    theButterflies->copyCurrentState(maxButterflyProfile);
    double maxButterfliesDensity = theButterflies->totalButterflyPopulation();
    double prevButterflyDensity = maxButterfliesDensity;
    double minButterfliesDensity = maxButterfliesDensity;
    int countButterflyIncreasing = 0;

    double prevWaspDensity = theButterflies->waspPopulation();
    double maxWaspDensity = prevWaspDensity;
    double minWaspDensity = maxWaspDensity;
    int countWaspIncreasing = 0;
    int prevValueClose = 0;


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
                if(fabs(maxButterfliesDensity-prevButterflyDensity)<1E-4)
                    prevValueClose = (prevValueClose | 0x1);
                else
                    prevValueClose = 0;

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
                    prevValueClose = (prevValueClose | 0x2);
                else
                    prevValueClose = 0;

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
                if(fabs(maxWaspDensity-prevWaspDensity)<1E-4)
                    prevValueClose = (prevValueClose | 0x4);
                else
                    prevValueClose = 0;

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
                    prevValueClose = (prevValueClose | 0x8);
                else
                    prevValueClose = 0;

                minWaspDensity = prevWaspDensity;
                countWaspIncreasing = 0;
            }
        }
        prevWaspDensity = currentWaspDensity;

        if((prevValueClose>=0xf)||(maxButterfliesDensity>10.0))
            break;
    }

    // Clean up the allocated space
    delete theButterflies;
    arrays.delonetensor(maxButterflyProfile);

    MaxMinBuffer values;
    values.which = which;
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
    msgsnd(msgID,&values,sizeof(values),0);
    //std::cout << "DONE " << which << std::endl;

    return(1);

}

int NumericalTrials::approximateSystemHysteresis(
        double mu,
        double c, double g, double d,
        double mLow, double mHigh, int numberM,
        double dt, int maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint,
        bool appendFile)
{

    // Open a file to write the results to.  Write the header for the file as well.
    std::fstream csvFile;

    if(appendFile)
    {
        csvFile.open("changingMHysteresisReverse.csv", std::ios::out | std::ios::app);
    }
    else
    {
        csvFile.open("changingMHysteresisReverse.csv", std::ios::out);
        csvFile << "mu,c,g,d,m,time,maxWasp,minWasp,minButterfly,maxButterfly" << std::endl;
    }

    Butterflies *theButterflies = new Butterflies(legendrePolyDegree,legendrePolyDegree+2);
    theButterflies->initializeLegendreParams();
    theButterflies->setMu(mu);
    theButterflies->setC(c);
    theButterflies->setG(g);
    theButterflies->setD(d);
    theButterflies->setDT(dt);

    theButterflies->initializeButterflies();
    //theButterflies->initializeButterfliesGaussian(1.0,mu*0.5);


    double maxButterfliesDensity = 0.0;
    double minButterfliesDensity = 0.0;
    double maxWaspDensity = 0.0;
    double minWaspDensity = 0.0;
    double timeSpan = 0.0;

    long mLupe;
    double currentM = mLow;
    for(mLupe=0;mLupe<=static_cast<long>(numberM);++mLupe)
    {
        currentM = mLow + static_cast<double>(mLupe)*(mHigh-mLow)/static_cast<double>(numberM-1);
        std::cout << "Starting approximation " << mLupe << ", " << currentM << std::endl;
        theButterflies->setM(currentM);
        approximateSystemGivenInitial(theButterflies,timeSpan,dt,maxTimeLupe,
                                      legendrePolyDegree,maxDeltaNorm,maxNewtonSteps,skipPrint,
                                      maxButterfliesDensity,minButterfliesDensity,
                                      maxWaspDensity,minWaspDensity);


        csvFile << theButterflies->getMu() << "," << theButterflies->getC() << "," << theButterflies->getG() << ","
                << theButterflies->getD() << "," << theButterflies->getM() << ","
                << timeSpan << ","
                << maxWaspDensity << "," << minWaspDensity << ","
                << minButterfliesDensity << "," << maxButterfliesDensity << std::endl;

    }
    std::cout << std::endl << std::endl << "All processes finished." << std::endl;


    delete theButterflies;
    csvFile.close();
    return(1);
}

int NumericalTrials::approximateSystemGivenInitial(Butterflies *theButterflies,
        double &timeSpan,double dt, int maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint,
        double &maxButterfliesDensity, double &minButterfliesDensity,
        double &maxWaspDensity, double &minWaspDensity)
{
    ArrayUtils<double> arrays;
    double *maxButterflyProfile = arrays.onetensor(legendrePolyDegree+2);
    int    timeLupe = 0;
    timeSpan        = 0.0;


    //std::cout << "Pre-processing" << std::endl;
    theButterflies->copyCurrentState(maxButterflyProfile);
    maxButterfliesDensity = theButterflies->totalButterflyPopulation();
    double prevButterflyDensity = maxButterfliesDensity;
    minButterfliesDensity = maxButterfliesDensity;
    int countButterflyIncreasing = 0;

    double prevWaspDensity = theButterflies->waspPopulation();
    maxWaspDensity = prevWaspDensity;
    minWaspDensity = maxWaspDensity;
    int countWaspIncreasing = 0;
    int prevValueClose = 0;


    // Start the time loop, and calculation an approximation at
    // each time step.
    for(timeLupe=0;timeLupe<maxTimeLupe;++timeLupe)
    {
        timeSpan = static_cast<double>(timeLupe)*dt;

        if((skipPrint>0) && (timeLupe%(skipPrint)==0))
        {
            std::cout << "Calculating an approximation: "
                         << std::fixed
                         << std::setw(8)
                         << std::setprecision(4)
                         << timeLupe << " (" << timeSpan <<  ") ";
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
                if(fabs(maxButterfliesDensity-prevButterflyDensity)<1E-4)
                    prevValueClose = (prevValueClose | 0x1);
                else
                    prevValueClose = 0;

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
                    prevValueClose = (prevValueClose | 0x2);
                else
                    prevValueClose = 0;

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
                if(fabs(maxWaspDensity-prevWaspDensity)<1E-4)
                    prevValueClose = (prevValueClose | 0x4);
                else
                    prevValueClose = 0;

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
                    prevValueClose = (prevValueClose | 0x8);
                else
                    prevValueClose = 0;

                minWaspDensity = prevWaspDensity;
                countWaspIncreasing = 0;
            }
        }
        prevWaspDensity = currentWaspDensity;

        if((prevValueClose>=0xf)||(maxButterfliesDensity>10.0))
            break;
    }

    arrays.delonetensor(maxButterflyProfile);

    return(1);

}

