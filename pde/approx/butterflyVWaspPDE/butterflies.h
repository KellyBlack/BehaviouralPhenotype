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
    virtual void calculateRHS();
    virtual void copyCurrentStateToTemp();
    virtual void copyCurrentState(double *ptr);
    double totalButterflyPopulation();
    double waspPopulation();

    void initializeButterflies();
    void initializeButterflies(double *initial);
    void initializeButterfliesGaussian(double center, double variance);
    void initializeButterfliesConstant(double theta);
    void deleteButterflies();
    virtual void writeParameters(std::fstream &resultsFile);
    virtual void writeCurrentApprox(double time,std::ofstream &resultsFile);
    virtual void writeBinaryCurrentApprox(double &time, std::fstream &resultsFile);

    void setMu(double value) { mu = value; }
    double getMu() { return(mu); }

    void setC(double value) { c = value; }
    double getC() { return(c); }

    void setG(double value) { g = value; }
    double getG() { return(g); }

    void setD(double value) { d = value; }
    double getD() { return(d); }

    void setM(double value) { m = value; }
    double getM() { return(m); }

    double getLeftEndButterflies();
    double getRightEndButterflies();
    double getButterfly(int which);
    void   setButterfly(int which,double val);
    void   copyState(Butterflies* butterfly);
    void   printState();

    void setLeftMax (double value) { leftMax  = value;}
    void setLeftMin (double value) { leftMin  = value;}
    void setRightMax(double value) { rightMax = value;}
    void setRightMin(double value) { rightMin = value;}

    double getLeftMax()  { return(leftMax);}
    double getLeftMin()  { return(leftMin);}
    double getRightMax() { return(rightMax);}
    double getRightMin() { return(rightMin);}

    void   setLastTime(double value) { lastTime = value;}
    double getLastTime() { return(lastTime);}

protected:

    double parameterDistribution(double theta);

private:

    double *butterflies  = nullptr;
    double *prevTimeStep = nullptr;

    double mu = 1.0;
    double c  = 3.0;
    double g  = 0.2;
    double d  = 0.5;
    double m  = 0.7;

    double leftMax  = 0.0;
    double leftMin  = 0.0;
    double rightMax = 0.0;
    double rightMin = 0.0;
    double lastTime = 0.0;
};

#endif // BUTTERFLIES_H
