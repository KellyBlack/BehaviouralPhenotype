#ifndef PDESOLVER_H
#define PDESOLVER_H

#include <fstream>

#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"


class PDESolver
{
public:
    PDESolver(int number=0);

    void setNumber(int number);
    int getNumber();

    void createArrays();
    void deleteArrays();
    void initializeLegendreParams();

    virtual void buildJacobian() = 0;
    virtual void updateNewtonStep() = 0;
    bool solveLinearizedSystem();
    double normDelta();

    void writeAbscissa(std::ofstream &resultsFile);
    void writeCurrentApprox(std::ofstream &resultsFile);

protected:

    ~PDESolver();

    int N;

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
    int    *order     = nullptr;

private:

};

#endif // PDESOLVER_H
