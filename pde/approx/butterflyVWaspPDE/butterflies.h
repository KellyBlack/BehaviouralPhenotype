#ifndef BUTTERFLIES_H
#define BUTTERFLIES_H

#include <fstream>

#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"


class butterflies
{
public:
    butterflies(int number=0);
    ~butterflies();

    void setNumber(int number);
    int getNumber();

    void createArrays();
    void deleteArrays();
    void initializeLegendreParams();

    void buildJacobian();
    bool solveLinearizedSystem();

    void writeAbscissa(std::ofstream &resultsFile);
    void writeCurrentApprox(std::ofstream &resultsFile);

protected:

private:
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

};

#endif // BUTTERFLIES_H
