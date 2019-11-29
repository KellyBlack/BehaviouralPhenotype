#ifndef NUMERICALTRIALS_H
#define NUMERICALTRIALS_H

#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"
#include "butterflies.h"


class NumericalTrials
{
public:
    NumericalTrials();

    int approximateSystem(
            double mu, double c, double g, double d, double m,
            double dt, int maxTimeLupe,
            int legendrePolyDegree,
            double maxDeltaNorm, int maxNewtonSteps,
            std::fstream &binFile,
            int skipPrint,
            int skipFileSave);

};

#endif // NUMERICALTRIALS_H
