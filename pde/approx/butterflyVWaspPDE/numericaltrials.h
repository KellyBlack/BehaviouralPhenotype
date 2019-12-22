#ifndef NUMERICALTRIALS_H
#define NUMERICALTRIALS_H

#include <vector>
#include <string>
#include <thread>


#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"
#include "butterflies.h"


class NumericalTrials
{
public:
    NumericalTrials();

    static void multipleApproximationsByM(
            double mu, double c, double g, double d,
            double lowM, double highM, double stepM,
            double dt, int maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint,
            int skipFileSave,
            int numberThreads);

    int approximateSystem(
            double mu, double c, double g, double d, double m,
            double dt, int maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            std::string filename,
            int skipPrint,
            int skipFileSave);

    int approximateSystemTrackRepeating(double muLow, double muHigh, int numberMu,
            double c, double g, double d,
            double mLow, double mHigh, int numberM,
            double dt, int maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint, int numProcesses, bool appendFile);

    int approximateSystemQuietResponse(double mu, double c, double g, double d, double m,
            double dt, int maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint, int msgID,unsigned long which);

    int approximateSystemHysteresis(double mu,
            double c, double g, double d,
            double mLow, double mHigh, int numberM,
            double dt, int maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint,
            bool appendFile);

    int approximateSystemGivenInitial(Butterflies *theButterflies,
            double &timeSpan, double dt, int maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint,
            double &maxButterfly, double &minButterfly,
            double &maxWasp, double &minWaspDensity);

protected:

    // Set up the data structure that will be used to keep track of the processes
    // and returned information associated with each process.
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

    struct MaxMinBuffer
    {
      long mtype;
      unsigned long which;
      double mu;
      double c;
      double g;
      double d;
      double m;
      double endTime;
      double maxButterfly;
      double minButterfly;
      double maxWasp;
      double minWasp;
    };

    MessageInformation* findReturnedProcessParameters(
            MaxMinBuffer msgValue,
            std::vector<MessageInformation*> processes);

private:



};

#endif // NUMERICALTRIALS_H
