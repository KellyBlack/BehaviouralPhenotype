#ifndef NUMERICALTRIALS_H
#define NUMERICALTRIALS_H

#include <string>

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

};

#endif // NUMERICALTRIALS_H
