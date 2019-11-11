#ifndef BUTTERFLIES_H
#define BUTTERFLIES_H

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

protected:

private:
    int N;

    double **lhsArray    = nullptr;
    double **secondDeriv = nullptr;
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
