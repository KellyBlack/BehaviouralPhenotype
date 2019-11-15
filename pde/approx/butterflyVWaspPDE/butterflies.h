#ifndef BUTTERFLIES_H
#define BUTTERFLIES_H

#include "pdeSolver.h"

class Butterflies : public PDESolver
{
public:
    Butterflies(int number=0,int sizeState=0);
    virtual ~Butterflies();

    virtual void buildJacobian();
    virtual void updateNewtonStep();
    void calculateRHS();
    void copyCurrentStateToTemp();

    void initializeButterflies();
    void deleteButterflies();
    virtual void writeCurrentApprox(double time,std::ofstream &resultsFile);

    void setMu(double value) { mu = value; }
    double getMu() { return(mu); }

    void setC(double value) { c = value; }
    double getC() { return(c); }

    void setF(double value) { f = value; }
    double getF() { return(f); }

    void setG(double value) { g = value; }
    double getG() { return(g); }

    void setD(double value) { d = value; }
    double getD() { return(d); }

protected:

    double parameterDistribution(double theta);

private:

    double *butterflies  = nullptr;
    double *prevTimeStep = nullptr;

    double mu = 1.0;
    double c  = 3.0;
    double f  = 2.0;
    double g  = 0.2;
    double d  = 0.5;
};

#endif // BUTTERFLIES_H
