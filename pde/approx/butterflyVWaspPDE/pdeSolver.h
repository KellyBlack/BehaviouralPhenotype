#ifndef PDESOLVER_H
#define PDESOLVER_H

#include <fstream>

#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"


class PDESolver
{
public:
    PDESolver(int number=0,int numberState=0);

    void setNumber(int number);
    int getNumber();

    void setStateSize(int number);
    int getStateSize();

    void setDT(double value);
    double getDT();

    void createArrays();
    void deleteArrays();
    void initializeLegendreParams();
    
    virtual void buildJacobianTimeStepping() = 0;
    virtual void buildjacobianSteadyState()  = 0;
    virtual void updateNewtonStep() = 0;
    virtual void calculateRHSTimeStepping() = 0;
    virtual void calculateRHSSteadyState() = 0;
    virtual void copyCurrentStateToTemp() = 0;
    virtual void copyCurrentState(double *ptr) = 0;
    void copyCurrentState(PDESolver *ptr);
    virtual void setState(int n,double value) = 0;
    virtual double getState(int n) = 0;
    bool solveLinearizedSystem();
    double normDelta();
    int singleTimeStep(double maxNewtonDiffNorm,int maxNewtonSteps,bool printInfo);
    int steadyStateApprox(double maxNewtonDiffNorm,int maxNewtonSteps,bool printInfo);

    void writeAbscissa(std::ofstream &resultsFile);
    void writeBinaryHeader(std::fstream &resultsFile);
    virtual void writeParameters(std::fstream &resultsFile) = 0;
    virtual void writeCurrentApprox(double time,std::ofstream &resultsFile) = 0;
    virtual void writeBinaryCurrentApprox(double &time,std::fstream &resultsFile) = 0;

protected:

    ~PDESolver();

    int N;
    int stateSize;
    double dt = 0.0;

    // Define the matrices used for
    double **lval        = nullptr;
    double **D1          = nullptr;
    double **stiff       = nullptr;

    // Define the vectors that are used.
    // Includes the Gauss quadrature (abscissa and weights)
    double *gaussWeights  = nullptr;
    double *gaussAbscissa = nullptr;

    // Define the matrices and vectors used to solve the linear systems.
    double **jacobian = nullptr;
    double *baseFunc  = nullptr;
    double *deltaX    = nullptr;
    double *rhs       = nullptr;
    int    *order     = nullptr;

private:

};

#endif // PDESOLVER_H
