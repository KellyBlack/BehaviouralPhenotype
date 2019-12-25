#ifndef RUNGAKUTTA45_H
#define RUNGAKUTTA45_H

#include <string>
#include <thread>

#include "approximationbase.h"

class RungaKutta45 : public ApproximationBase
{
public:
    explicit RungaKutta45(double initialTimeStep=1.0E-3,double *initialCondition=nullptr);


    long approximationByM(double cValue, double gValue, double dValue, double thetaValue,
                      double lowM, double highM, long numberM,
                      double startTime, double endTime, double initialDt, double minimumDT, double *initialCond, double tolerance,
                      std::string filename, bool appendFile,
                      int numberThreads);

    long approximation(long which,
                      double cValue, double gValue, double dValue, double mValue, double thetaValue,
                      double startTime, double endTime, double initialDt, double minimumDT,
                      double *initialCond,double tolerance, int msgID);


    void setC(double cValue) { c = cValue; }
    double getC() { return(c); }

    void setG(double gValue) { g = gValue; }
    double getG() { return(g); }

    void setD(double dValue) { d = dValue; }
    double getD() { return(d); }

    void setM(double mValue) { m = mValue; }
    double getM() { return(m); }

    void setTheta(double thetaValue) { theta = thetaValue; }
    double getTheta() { return(theta); }

    void setDT(double value) { dt = value; }
    double getDT() { return(dt); }

    void setCurrentTime(double value) { currentTime = value; }
    double getCurrentTime() { return(currentTime); }

    void setMinimumDT(double value) { minimumDT = value; }
    double getMinimumDT() { return(minimumDT); }

protected:

    void rateFunction(double theTime,double *x,double *rate);
    double singleStep(double tolerance);

private:
    double magnitudeDifference(double *one,double *two);
    double state[2];
    double dt;
    double currentTime;

    double c = 0.0;
    double g = 0.0;
    double d = 0.0;
    double m = 0.0;
    double theta = 0.0;
    double minimumDT = 1.0E-5;



};

#endif // RUNGAKUTTA45_H
