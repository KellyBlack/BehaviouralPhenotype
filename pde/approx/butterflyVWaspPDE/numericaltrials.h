#ifndef NUMERICALTRIALS_H
#define NUMERICALTRIALS_H

#include <vector>
#include <string>


#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"
#include "approximationbase.h"
#include "butterflies.h"


class NumericalTrials : public ApproximationBase
{
public:
    NumericalTrials();

    static void multipleApproximationsByM(double mu, double c, double g, double d,
            double lowM, double highM, double stepM,
            double dt, unsigned long maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint,
            int skipFileSave,
            int numberThreads);

    static void multipleApproximationsByMandC(double mu, double g, double d,
            double lowC, double highC, double stepC,
            double lowM, double highM, double stepM,
            double dt, unsigned long maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint, int skipFileSave,
            int numberThreads,
            std::string filename);

    int approximateSystem(
            double mu, double c, double g, double d, double m,
            double dt,unsigned long maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            std::string filename,
            int skipPrint,
            int skipFileSave);

    int approximateSystemCheckOscillation(
            double mu, double c, double g, double d, double m,
            double dt, unsigned long maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            std::string filename,
            int skipPrint
            );

    int approximateSystemTrackRepeating(double muLow, double muHigh, int numberMu,
            double c, double g, double d,
            double mLow, double mHigh, int numberM,
            double dt, unsigned long maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint, int numProcesses, bool appendFile, std::string filename);

    int approximateSystemQuietResponse(
            double mu, double c, double g, double d, double m,
            double dt, unsigned long maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint, int msgID, long which);

    int approximateSystemHysteresis(double mu,
            double c, double g, double d,
            double mLow, double mHigh, int numberM,
            double dt,unsigned long maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint,
            bool appendFile);

    int approximateSystemGivenInitial(Butterflies *theButterflies,
            double &timeSpan, double dt,unsigned long maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            int skipPrint,
            double &maxButterfly, double &minButterfly,
            double &maxWasp, double &minWaspDensity);

protected:


private:

    bool checkRepeating(
            Butterflies *theButterflies,
            double &prevButterflyDensity,
            double &maxButterfliesDensity,
            double &minButterfliesDensity,
            double &prevButterflyCheck,
            double &prevWaspDensity,
            double &maxWaspDensity,
            double &minWaspDensity,
            double &prevWaspCheck,
            int &prevCycleClose,
            int &prevValueClose,
            int &countButterflyIncreasing,
            int &countWaspIncreasing
            );



};

#endif // NUMERICALTRIALS_H
