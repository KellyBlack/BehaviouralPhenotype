#include <fstream>
#include <string>
#include <sstream>
#include <unistd.h>
#include <atomic>

#include <sys/ipc.h>
#include <sys/msg.h>


#include "numericaltrials.h"
#include "util.h"
#include "limInf.h"

//#define OUTPUTFILE "approximation.csv"
#define BINARYOUTPUTFILE "approximation-TRIAL-"


// Constructor for the numerical trials class.
NumericalTrials::NumericalTrials() :
    ApproximationBase::ApproximationBase()
{
    theButterflies = nullptr;
}

void NumericalTrials::multipleApproximationsByM(
        double mu, double c, double g, double d,
        double lowM, double highM, double stepM,
        double dt, unsigned long maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint, int skipFileSave,
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
                     << "-c-" <<  std::setw(6) << std::fixed << std::setprecision(4) << std::setfill('0') << c
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

void NumericalTrials::multipleApproximationsByMandC(
        double mu,double g, double d,
        double lowC, double highC, double stepC,
        double lowM, double highM, double stepM,
        double dt, unsigned long maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint, int skipFileSave,
        int numberThreads,
        std::string filename)
{
    // Open a file to write the results to.  Write the header for the file as well if this is a new file..
    std::fstream csvFile;
    csvFile.open(filename, std::ios::out);
    csvFile << "which,mu,c,g,d,m,time,maxButterflyLeft,minButterflyLeft,maxButterflyRight,minButterflyRight" << std::endl;
    csvFile.close();

    double m = lowM;
    double c = lowC;

    NumericalTrials* lastTrial = nullptr;

    struct remoteProcess {
        std::thread process;
        std::atomic<bool> running;
        NumericalTrials* trial;
    };

    std::vector<remoteProcess*> processes;

    processes.clear();
    while((m<=highM) || (!processes.empty()))
    {
        while((m<=highM)&&(static_cast<int>(processes.size())<numberThreads))
        {
            NumericalTrials *trial    = new NumericalTrials();
            remoteProcess* newProcess = new remoteProcess;
            newProcess->running = true;
            newProcess->trial = trial;
            newProcess->process = std::thread(&NumericalTrials::approximateSystemCheckOscillation,trial,
                                              mu,c,g,d,m,
                                              dt,maxTimeLupe,
                                              legendrePolyDegree,
                                              maxDeltaNorm,maxNewtonSteps,
                                              filename,
                                              skipPrint,skipFileSave,lastTrial,&(newProcess->running));

            processes.push_back(newProcess);
            c += stepC;
            if(c>highC)
            {
                c = lowC;
                m += stepM;
            }

        }

        sleep(2);
        std::cout << "Checking processes to see if any have stopped, " << processes.size() << "." << std::endl;

        std::vector<remoteProcess*>::iterator eachProcess;
        for(eachProcess=processes.begin();eachProcess!=processes.end();++eachProcess)
        {

            remoteProcess* process = *eachProcess;
            if(!process->running)
            {
                std::cout << "Thread has finished: " << process << std::endl;
                if(process->process.joinable())
                {
                    if(lastTrial!=nullptr)
                    {
                        delete lastTrial;
                    }
                    //lastTrial = process->trial;
                    delete process->trial;

                    process->process.join();
                    process->trial = nullptr;
                    processes.erase(eachProcess);
                    delete process;

                    //eachProcess--;
                    break;
                }
            }

        }

    }

    processes.clear();

}

int NumericalTrials::approximateSystem(
        double mu, double c, double g, double d, double m,
        double dt, unsigned long maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        std::string filename,
        int skipPrint, int skipFileSave)
{
    double t        = 0.0;
    unsigned long timeLupe = 0;

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

    //theButterflies->initializeButterfliesGaussian(1.0,mu);
    theButterflies->initializeButterfliesConstant(1.0);
    theButterflies->writeParameters(binFile);
    theButterflies->writeBinaryHeader(binFile);

    // Define the variables used to determing if the system is repeating
    // Used to figure out when to stop.
    double maxButterfliesDensity     = theButterflies->totalButterflyPopulation();
    double prevButterflyDensity      = maxButterfliesDensity;
    double minButterfliesDensity     = maxButterfliesDensity;
    double prevMaxButterfliesDensity = maxButterfliesDensity;
    double prevMinButterfliesDensity = maxButterfliesDensity;
    int countButterflyIncreasing = 0;

    double prevWaspDensity    = theButterflies->waspPopulation();
    double maxWaspDensity     = prevWaspDensity;
    double minWaspDensity     = maxWaspDensity;
    double prevMaxWaspDensity = maxWaspDensity;
    double prevMinWaspDensity = maxWaspDensity;
    int countWaspIncreasing   = -2;
    int prevCycleClose = 0;

    double prevButterflyCheck = 0.0;
    double prevWaspCheck = 0.0;
    int prevValueClose = 0;

    LimInf<double> maxButterflyLeft(theButterflies->getLeftThirdButterflies(),true);
    LimInf<double> minButterflyLeft(theButterflies->getLeftThirdButterflies(),false);
    LimInf<double> maxButterflyRight(theButterflies->getRightThirdButterflies(),true);
    LimInf<double> minButterflyRight(theButterflies->getRightThirdButterflies(),false);


    // Start the time loop, and calculation an approximation at
    // each time step.
    for(timeLupe=0;timeLupe<maxTimeLupe;++timeLupe)
    {
        t = static_cast<double>(timeLupe)*dt;

        if(timeLupe%(static_cast<unsigned long>(skipPrint))==0)
        {
            std::cout << "Calculating an approximation: "
                         << std::fixed
                         << std::setw(8)
                         << std::setprecision(4)
                         << timeLupe << " (" << t << ") "
                         << maxButterflyLeft.extreme() << "  "
                         << minButterflyLeft.extreme() << "  "
                         << maxButterflyRight.extreme() << "  "
                         << minButterflyRight.extreme() << "  "
                            ;
        }

        if(
           theButterflies->singleTimeStep(maxDeltaNorm,maxNewtonSteps,timeLupe%(static_cast<unsigned long>(skipPrint))==0)
                < 0)
        {
            std::cout << std::endl << "Error - Newton's Method did not converge." << std::endl;
            return(0);
        }

        if(timeLupe%static_cast<unsigned long>(skipFileSave)==0)
        {
            theButterflies->writeBinaryCurrentApprox(t,binFile);
        }

        int repeating = checkRepeating(
                    theButterflies,
                    prevButterflyDensity,maxButterfliesDensity,minButterfliesDensity,prevMaxButterfliesDensity,prevMinButterfliesDensity,prevButterflyCheck,
                    prevWaspDensity,maxWaspDensity,minWaspDensity,prevMaxWaspDensity,prevMinWaspDensity,prevWaspCheck,
                    prevCycleClose,prevValueClose,countButterflyIncreasing,countWaspIncreasing
                    );
        if(repeating != 0)
        {
            if(repeating > 0)
            {
                maxButterflyLeft.setExtreme(theButterflies->getLeftThirdButterflies());
                minButterflyLeft.setExtreme(theButterflies->getLeftThirdButterflies());
                maxButterflyRight.setExtreme(theButterflies->getRightThirdButterflies());
                minButterflyRight.setExtreme(theButterflies->getRightThirdButterflies());
                std::cout << "Steady state achieved." << std::endl;
            }
            break; // the solution is either settling into a steady state or repeating.... enough is enough.
        }

    }

    // Clean up the data file and close it
    binFile.close();
    delete theButterflies;

    return(1);
}

int NumericalTrials::approximateSystemCheckOscillation(double mu, double c, double g, double d, double m,
        double dt, unsigned long maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        std::string filename,
        int skipPrint, int skipFileSave,
        NumericalTrials *prevApprox,
        std::atomic<bool> *running
        )
{

    *running               = true;
    double t               = 0.0;
    unsigned long timeLupe = 0;

    // Variables used to save the results of calculations into a
    // data file.

    std::cout << "Pre-processing" << std::endl;
    int N = legendrePolyDegree;
    if(theButterflies==nullptr)
        theButterflies = new Butterflies(N,N+2);
    theButterflies->initializeLegendreParams();
    theButterflies->setMu(mu);
    theButterflies->setC(c);
    theButterflies->setG(g);
    theButterflies->setD(d);
    theButterflies->setM(m);
    theButterflies->setDT(dt);

    //theButterflies->initializeButterfliesGaussian(1.0,mu);
    theButterflies->initializeButterfliesConstant(1.0);
    //theButterflies->writeParameters(binFile);
    //theButterflies->writeBinaryHeader(binFile);

    if(prevApprox != nullptr)
        theButterflies->copyState(prevApprox->getButterflies());

    // Define the variables used to determing if the system is repeating
    // Used to figure out when to stop.
    double maxButterfliesDensity = theButterflies->totalButterflyPopulation();
    double prevButterflyDensity = maxButterfliesDensity;
    double minButterfliesDensity = maxButterfliesDensity;
    double prevMinButterfliesDensity = maxButterfliesDensity;
    double prevMaxButterfliesDensity = maxButterfliesDensity;
    int countButterflyIncreasing = 0;

    double prevWaspDensity    = theButterflies->waspPopulation();
    double maxWaspDensity     = prevWaspDensity;
    double minWaspDensity     = maxWaspDensity;
    double prevMaxWaspDensity = prevWaspDensity;
    double prevMinWaspDensity = maxWaspDensity;
    int countWaspIncreasing = -2;
    int prevCycleClose = 0;

    double prevButterflyCheck = 0.0;
    double prevWaspCheck = 0.0;
    int prevValueClose = 0;

    LimInf<double> maxButterflyLeft(theButterflies->getLeftThirdButterflies(),true);
    LimInf<double> minButterflyLeft(theButterflies->getLeftThirdButterflies(),false);
    LimInf<double> maxButterflyRight(theButterflies->getRightThirdButterflies(),true);
    LimInf<double> minButterflyRight(theButterflies->getRightThirdButterflies(),false);


    // Start the time loop, and calculation an approximation at
    // each time step.
    for(timeLupe=0;timeLupe<maxTimeLupe;++timeLupe)
    {
        t = static_cast<double>(timeLupe)*dt;

        if(timeLupe%(static_cast<unsigned long>(skipPrint))==0)
        {
            std::cout << "Calculating an approximation: "
                         << std::fixed
                         << std::setw(8)
                         << std::setprecision(4)
                         << timeLupe << " (" << t << ") "
                         << maxButterflyLeft.extreme() << "  "
                         << minButterflyLeft.extreme() << "  "
                         << maxButterflyRight.extreme() << "  "
                         << minButterflyRight.extreme() << "  "
                         ;
        }

        if(
           theButterflies->singleTimeStep(maxDeltaNorm,maxNewtonSteps,timeLupe%(static_cast<unsigned long>(skipPrint))==0)
                < 0)
        {
            std::cout << std::endl << "Error - Newton's Method did not converge." << std::endl;
            *running = false;
            return(0);
        }

        int repeating = checkRepeating(
                    theButterflies,
                    prevButterflyDensity,maxButterfliesDensity,minButterfliesDensity,prevMaxButterfliesDensity,prevMinButterfliesDensity,prevButterflyCheck,
                    prevWaspDensity,maxWaspDensity,minWaspDensity,prevMaxWaspDensity,prevMinWaspDensity,prevWaspCheck,
                    prevCycleClose,prevValueClose,countButterflyIncreasing,countWaspIncreasing
                    );
        if(repeating != 0)
        {
            if(repeating > 0)
            {
                maxButterflyLeft.setExtreme(theButterflies->getLeftThirdButterflies());
                minButterflyLeft.setExtreme(theButterflies->getLeftThirdButterflies());
                maxButterflyRight.setExtreme(theButterflies->getRightThirdButterflies());
                minButterflyRight.setExtreme(theButterflies->getRightThirdButterflies());
                maxWaspDensity = theButterflies->waspPopulation();
                minWaspDensity = theButterflies->waspPopulation();
            }
            break; // the solution is either settling into a steady state or repeating.... enough is enough.
        }

    }

    // Open a file to write the results to.  Write the header for the file as well if this is a new file..
    if(skipFileSave!=0)
    {
        std::fstream csvFile;
        csvFile.open(filename, std::ios::out | std::ios::app);

        csvFile << "which," << mu << ","
                << c << ","
                << g << ","
                << d << ","
                << m << ","
                << t << ","
                << maxButterflyLeft.extreme() << ","
                << minButterflyLeft.extreme() << ","
                << maxButterflyRight.extreme() << ","
                << minButterflyRight.extreme() << std::endl;

        // Clean up the data file and close it
        csvFile.close();
    }

    delete theButterflies;
    *running = false;
    return(1);
}

int NumericalTrials::checkRepeating(
        Butterflies *theButterflies,
        double &prevButterflyDensity,
        double &maxButterfliesDensity,
        double &minButterfliesDensity,
        double &prevMaxButterfliesDensity,
        double &prevMinButterfliesDensity,
        double &prevButterflyCheck,
        double &prevWaspDensity,
        double &maxWaspDensity,
        double &minWaspDensity,
        double &prevMaxWaspDensity,
        double &prevMinWaspDensity,
        double &prevWaspCheck,
        int &prevCycleClose,
        int &prevValueClose,
        int &countButterflyIncreasing,
        int &countWaspIncreasing
        )
{
    const double CLOSE_STEADY_STATE = 1.0E-6;
    const double CLOSE_SMALL_CHANGE = 1.0E-4;
    const long MAX_CHECK_CONSTANT   = 2000;
    double currentButterflyDensity  = theButterflies->totalButterflyPopulation();
    double currentWaspDensity       = theButterflies->waspPopulation();

    // First check to see if the solution is settling into a stable steady state.
    if(((fabs(currentWaspDensity-prevWaspCheck)<CLOSE_STEADY_STATE)&&(fabs(currentButterflyDensity-prevButterflyCheck)<CLOSE_STEADY_STATE)))
    {
        // The last approximation is close to the last check state.
        prevValueClose += 1;
        //if(prevValueClose>1900)
        //    std::cout << "butterflies and wasps close " << prevValueClose << std::endl;
    }
     else
    {
        // The last approximation is away from the check state. Reset the check state and the counter.
        prevValueClose = 0;
        prevWaspCheck = currentWaspDensity;
        prevButterflyCheck = currentButterflyDensity;
    }


    if(currentButterflyDensity<prevButterflyDensity)
    {
        // The butterfly population is decreasing. If countButterflyIncreasing
        // is negative then this is the first time it has gone through this
        // check
        if(countButterflyIncreasing < 0)
        {
            //std::cout << which << ": max butterfly " << prevButterflyDensity << "-" << maxButterfliesDensity << std::endl;
            if(fabs(maxButterfliesDensity-prevMaxButterfliesDensity)<CLOSE_SMALL_CHANGE)
                prevCycleClose = (prevCycleClose | 0x1);
            else
                prevCycleClose = 0;

            prevMaxButterfliesDensity = maxButterfliesDensity;
            countButterflyIncreasing = 2;
        }
        else
        {
            countButterflyIncreasing = 1;
            minButterfliesDensity = currentButterflyDensity;
        }

    }
    else
    {
        // The butterfly population is increasing. If countButterflyIncreasing
        // is positive then this is the first time it has gone through this check.
        if(countButterflyIncreasing > 0)
        {
            //std::cout << which << ": min butterfly " << prevButterflyDensity << "-" << minButterfliesDensity << std::endl;
            if(fabs(minButterfliesDensity-prevMinButterfliesDensity)<CLOSE_SMALL_CHANGE)
                prevCycleClose = (prevCycleClose | 0x2);
            else
                prevCycleClose = 0;

            prevMinButterfliesDensity = minButterfliesDensity;
            countButterflyIncreasing = -2;
        }
        else
        {
            countButterflyIncreasing = -1;
            maxButterfliesDensity = prevButterflyDensity;
        }

    }
    prevButterflyDensity = currentButterflyDensity;


    if(currentWaspDensity<prevWaspDensity)
    {
        // The wasp population is decreasing. If countIncreasing is
        // pos. that means that the wasp count had a long trend of
        // increasing.
        if(countWaspIncreasing < 0)
        {
            //std::cout << which << ": max wasp " << prevWaspDensity << "-" << maxWaspDensity << std::endl;
            if(fabs(maxWaspDensity-prevMaxWaspDensity)<CLOSE_SMALL_CHANGE)
                prevCycleClose = (prevCycleClose | 0x4);
            else
                prevCycleClose = 0;

            prevMaxWaspDensity = maxWaspDensity;
            countWaspIncreasing = 2;
        }
        else
        {
            countWaspIncreasing = 1;
            minWaspDensity = currentWaspDensity;
        }

    }
    else
    {
        // The wasp population is increasing. If countIncreasing is
        // neg. that means that the wasp count had a long trend of
        // decreasing.
        if(countWaspIncreasing > 0)
        {
            //std::cout << which << ": min wasp " << prevWaspDensity << "-" << minWaspDensity << std::endl;
            if(fabs(minWaspDensity-prevMinWaspDensity)<CLOSE_SMALL_CHANGE)
                prevCycleClose = (prevCycleClose | 0x8);
            else
                prevCycleClose = 0;

            prevMinWaspDensity = minWaspDensity;
            countWaspIncreasing = -2;
        }
        else
        {
            countWaspIncreasing = -1;
            maxWaspDensity = prevWaspDensity;
        }
    }
    prevWaspDensity = currentWaspDensity;

    if((prevCycleClose>=0xf)||(maxButterfliesDensity>10.0))
    {
        //std::cout << "cycle! " << prevCycleClose << std::endl;
        return(-1);
    }
    else if (prevValueClose>MAX_CHECK_CONSTANT)
    {
        // The solution is assumed to be close to a steady state.
        //std::cout << "Steady state! " << currentButterflyDensity << " " << currentWaspDensity << std::endl;
        maxButterfliesDensity = currentButterflyDensity;
        minButterfliesDensity = currentButterflyDensity;
        maxWaspDensity = currentWaspDensity;
        minWaspDensity = currentWaspDensity;

        prevMaxButterfliesDensity = currentButterflyDensity;
        prevMaxButterfliesDensity = currentButterflyDensity;
        prevMaxWaspDensity = currentWaspDensity;
        prevMinWaspDensity = currentWaspDensity;
        return(1);
    }

    return(0);
}


int NumericalTrials::approximateSystemTrackRepeating(
        double muLow, double muHigh, int numberMu,
        double c, double g, double d,
        double mLow, double mHigh, int numberM,
        double dt, unsigned long maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint, int numProcesses,
        bool appendFile,
        std::string filename)
{

    // Create an ID and then create a message queue that will be associated with the ID
    key_t msg_key = ftok("butterfly",1995);
    int msgID = msgget(msg_key,0666|IPC_CREAT);

    // Open a file to write the results to.  Write the header for the file as well if this is a new file..
    std::fstream csvFile;
    if(appendFile)
    {
        csvFile.open(filename, std::ios::out | std::ios::app);
    }
    else
    {
        csvFile.open(filename, std::ios::out);
        csvFile << "which,mu,c,g,d,m,time,maxWasp,minWasp,minButterfly,maxButterfly" << std::endl;
    }

    double currentM = mLow;
    double currentDiffusion = muLow;

    // Set up the vector used to keep track of all the information that the
    // remote processes create.
    MaxMinBuffer msgValue;
    NumericalTrials::MessageInformation* processInformation;
    std::vector<NumericalTrials::MessageInformation*> processes;

    // Need to make sure there are no messages in the buffer left
    // over from previous runs that were prematurely terminated.
    // (I will not go into the details about the agony of learning this lesson.)
    while(msgrcv(msgID,&msgValue,sizeof(msgValue),2,IPC_NOWAIT)>0)
    {
        std::cout << "Previous message in the queue: " << msgValue.which << std::endl;
    }

    // Loop through all possible values of the diffusion and value of m.
    unsigned long mLupe;
    unsigned long diffusionLupe;
    long currentNumberProcesses = 0;
    long totalRuns = 0;
    for(mLupe=0;mLupe<static_cast<unsigned long>(numberM);++mLupe)
        for(diffusionLupe=0;diffusionLupe<static_cast<unsigned long>(numberMu);++diffusionLupe)
        {
            // Calculate the value of m and the diffusion coefficient.
            currentM = mLow + static_cast<double>(mLupe)*(mHigh-mLow)/static_cast<double>(numberM-1);
            currentDiffusion = muLow + static_cast<double>(diffusionLupe)*(muHigh-muLow)/static_cast<unsigned long>(numberMu-1);

            // Set up the data structure that will keep track of the values of the coefficient for this numerical trial.
            // Then start a new thread.
            NumericalTrials *newTrial = new NumericalTrials();
            NumericalTrials::MessageInformation *newProcess = new NumericalTrials::MessageInformation;
            newProcess->which = totalRuns;
            newProcess->mu    = currentDiffusion;
            newProcess->c     = c;
            newProcess->g     = g;
            newProcess->d     = d;
            newProcess->m     = currentM;
            newProcess->trial = static_cast<ApproximationBase*>(newTrial);
            newProcess->process = new std::thread(
                            &NumericalTrials::approximateSystemQuietResponse,newTrial,
                            currentDiffusion,c,g,d,currentM,dt,maxTimeLupe,
                            legendrePolyDegree,maxDeltaNorm,maxNewtonSteps,skipPrint,msgID,totalRuns
                            );
            processes.push_back(newProcess);
            std::cout << "starting " << currentDiffusion << "/" << currentM << ": " << totalRuns << "  " << newProcess->process << std::endl;

            totalRuns += 1;
            if(++currentNumberProcesses>=static_cast<long>(numProcesses))
            {
                // There are too many processes running. Need to wait for one to stop
                // before starting a new process.
                msgrcv(msgID,&msgValue,sizeof(msgValue),2,0);

                // One just ended. Figure out the values that were sent and record the information from the run.
                std::cout << msgValue.which << "--"
                          << msgValue.mu << "," << msgValue.c << "," << msgValue.g << ","
                          << msgValue.d << "," << msgValue.m << "," << msgValue.endTime << ","
                          << msgValue.maxWasp << "," << msgValue.minWasp << ","
                          << msgValue.minButterfly << "," << msgValue.maxButterfly << std::endl;

                processInformation = findReturnedProcessParameters(msgValue,processes);
                if(processInformation!= nullptr)
                {
                    csvFile << processInformation->which << ","
                            << processInformation->mu << "," << processInformation->c << "," << processInformation->g << ","
                            << processInformation->d  << "," << processInformation->m << "," << processInformation->time << ","
                            << processInformation->maxWasp << "," << processInformation->minWasp << ","
                            << processInformation->minButterfly << "," << processInformation->maxButterfly << std::endl;

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
        msgrcv(msgID,&msgValue,sizeof(msgValue),2,0);
        std::cout << msgValue.which << " "
                  << msgValue.mu << "," << msgValue.c << "," << msgValue.g << ","
                  << msgValue.d << "," << msgValue.m << "," << msgValue.endTime << ","
                  << msgValue.maxWasp << "," << msgValue.minWasp << ","
                  << msgValue.minButterfly << "," << msgValue.maxButterfly << std::endl;

        processInformation = findReturnedProcessParameters(msgValue,processes);
        if(processInformation != nullptr)
        {
            csvFile << processInformation->which << ","
                    << processInformation->mu << "," << processInformation->c << "," << processInformation->g << ","
                    << processInformation->d << ","  << processInformation->m << "," << processInformation->time << ","
                    << processInformation->maxWasp << "," << processInformation->minWasp << ","
                    << processInformation->minButterfly << "," << processInformation->maxButterfly << std::endl;

            delete processInformation->trial;
            processInformation->process->join();
            currentNumberProcesses -= 1;
        }
    }

    std::vector<NumericalTrials::MessageInformation*>::iterator eachProcess;
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
        double dt, unsigned long maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint, int msgID,long which)
{
    const long MAX_CHECK_CONSTANT = 2000;
    ArrayUtils<double> arrays;
    double *maxButterflyProfile = arrays.onetensor(legendrePolyDegree+2);
    double t        = 0.0;
    unsigned long timeLupe = 0;


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

    //theButterflies->initializeButterflies();
    //theButterflies->initializeButterfliesGaussian(1.0,mu*0.25);
    theButterflies->initializeButterfliesConstant(1.0);
    theButterflies->copyCurrentState(maxButterflyProfile);
    double maxButterfliesDensity = theButterflies->totalButterflyPopulation();
    double prevButterflyDensity = maxButterfliesDensity;
    double minButterfliesDensity = maxButterfliesDensity;
    int countButterflyIncreasing = 0;

    double prevWaspDensity = theButterflies->waspPopulation();
    double maxWaspDensity = prevWaspDensity;
    double minWaspDensity = maxWaspDensity;
    int countWaspIncreasing = 0;
    int prevCycleClose = 0;

    double prevButterflyCheck = 0.0;
    double prevWaspCheck = 0.0;
    int prevValueClose = 0;


    // Start the time loop, and calculation an approximation at
    // each time step.
    for(timeLupe=0;timeLupe<maxTimeLupe;++timeLupe)
    {
        t = static_cast<double>(timeLupe)*dt;

        if((skipPrint>0) && (timeLupe%(static_cast<unsigned long>(skipPrint))==0))
        {
            std::cout << "Calculating an approximation: "
                         << std::fixed
                         << std::setw(8)
                         << std::setprecision(4)
                         << timeLupe << " (" << t << "-" << which <<  ") ";
        }

        if(
           theButterflies->singleTimeStep(maxDeltaNorm,maxNewtonSteps,(skipPrint>0)&&(timeLupe%(static_cast<unsigned long>(skipPrint))==0))
                < 0)
        {
            std::cout << std::endl << "Error - Newton's Method did not converge." << std::endl;
            return(0);
        }

        double currentButterflyDensity = theButterflies->totalButterflyPopulation();
        double currentWaspDensity = theButterflies->waspPopulation();

        if((fabs(currentWaspDensity-prevWaspDensity)<1.0E-7)&&(fabs(currentButterflyDensity-prevButterflyDensity)<1.0E-7))
        {
            if(((fabs(prevWaspCheck-currentWaspDensity)>1.0E-7)||(fabs(prevButterflyCheck-currentButterflyDensity)>1.0E-7)))
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
    msgsnd(msgID,&values,sizeof(values),2);
    //std::cout << "DONE " << which << std::endl;

    return(1);

}

int NumericalTrials::approximateSystemHysteresis(double mu,
        double c, double g, double d,
        double mLow, double mHigh, int numberM,
        double dt, unsigned long maxTimeLupe,
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

    //theButterflies->initializeButterflies();
    //theButterflies->initializeButterfliesGaussian(1.0,mu*0.5);
    theButterflies->initializeButterfliesConstant(1.0);

    double maxButterfliesDensity = 0.0;
    double minButterfliesDensity = 0.0;
    double maxWaspDensity = 0.0;
    double minWaspDensity = 0.0;
    double timeSpan = 0.0;

    long mLupe;
    double currentM = mLow;
    for(mLupe=static_cast<long>(numberM);mLupe>=0;--mLupe)
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
        double &timeSpan, double dt, unsigned long maxTimeLupe,
        int legendrePolyDegree,
        double maxDeltaNorm, int maxNewtonSteps,
        int skipPrint,
        double &maxButterfliesDensity, double &minButterfliesDensity,
        double &maxWaspDensity, double &minWaspDensity)
{
    const long MAX_CHECK_CONSTANT = 2000;
    ArrayUtils<double> arrays;
    double *maxButterflyProfile = arrays.onetensor(legendrePolyDegree+2);
    unsigned long timeLupe = 0;
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
    int prevCycleClose = 0;

    double prevButterflyCheck = 0.0;
    double prevWaspCheck = 0.0;
    int prevValueClose = 0;


    // Start the time loop, and calculation an approximation at
    // each time step.
    for(timeLupe=0;timeLupe<maxTimeLupe;++timeLupe)
    {
        timeSpan = static_cast<double>(timeLupe)*dt;

        if((skipPrint>0) && (timeLupe%(static_cast<unsigned long>(skipPrint))==0))
        {
            std::cout << "Calculating an approximation: "
                         << std::fixed
                         << std::setw(8)
                         << std::setprecision(4)
                         << timeLupe << " (" << timeSpan <<  ") ";
        }

        if(
           theButterflies->singleTimeStep(maxDeltaNorm,maxNewtonSteps,(skipPrint>0)&&(timeLupe%(static_cast<unsigned long>(skipPrint))==0))
                < 0)
        {
            std::cout << std::endl << "Error - Newton's Method did not converge." << std::endl;
            return(0);
        }

        double currentButterflyDensity = theButterflies->totalButterflyPopulation();
        double currentWaspDensity = theButterflies->waspPopulation();

        if((fabs(currentWaspDensity-prevWaspDensity)<1.0E-7)&&(fabs(currentButterflyDensity-prevButterflyDensity)<1.0E-7))
        {
            if(((fabs(prevWaspCheck-currentWaspDensity)>1.0E-7)||(fabs(prevButterflyCheck-currentButterflyDensity)>1.0E-7)))
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
                if((prevCycleClose>0x4)&&(fabs(maxButterfliesDensity-prevButterflyDensity)<1E-4))
                    prevCycleClose = (prevCycleClose | 0x8);
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
                    prevCycleClose = (prevCycleClose | 0x4);
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
                    prevCycleClose = (prevCycleClose | 0x2);
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
                    prevCycleClose = (prevCycleClose | 0x1);
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

    arrays.delonetensor(maxButterflyProfile);

    return(1);

}

