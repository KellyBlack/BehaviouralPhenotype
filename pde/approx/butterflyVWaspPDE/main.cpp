#include <iostream>
#include <fstream>
#include <cmath>

#include "util.h"
#include "legendre.h"
#include "lu_decomp.h"

#define OUTPUTFILE "approximation.csv"

int main()
{
    std::cout << "Pre-processing starting" << std::endl;

    //ArrayUtils<double> arrayCreator;
    //Legendre<double> legendre;

    // Define the number of grid points to use.
    // Define the matrices that are used to create the system of equations.
    int N = 40;
    double **lhsArray    = ArrayUtils<double>::twotensor(N+1,N+1);
    double **secondDeriv = ArrayUtils<double>::twotensor(N+1,N+1);
    double **lval        = ArrayUtils<double>::twotensor(N+1,N+1);
    double **D1          = ArrayUtils<double>::twotensor(N+1,N+1);
    double **stiff       = ArrayUtils<double>::twotensor(N+1,N+1);

    // Define the vectors that are used.
    // Includes the Gauss quadrature (abscissa and weights)
    double *gaussWeights  = ArrayUtils<double>::onetensor(N+1);
    double *gaussAbscissa = ArrayUtils<double>::onetensor(N+1);

    // Define the matrices and vectors used to solve the linear systems.
    double **jacobian = ArrayUtils<double>::twotensor(N+1,N+1);
    double *baseFunc  = ArrayUtils<double>::onetensor(N+1);
    double *deltaX  = ArrayUtils<double>::onetensor(N+1);
    int    *order     = ArrayUtils<int>::onetensor(N+1);

    // Miscellaneous helper variables, for loops for example.
    int innerLupe;
    int outerLupe;

    // Define the Gauss quadrature.
    Legendre<double>::leg_quad(gaussAbscissa,gaussWeights,N);

    // Define the stiffness and mass matrices for the Legendre collocation method.
    std::cout << "Pre-processing" << std::endl;
    Legendre<double>::leg_val(lval,gaussAbscissa,N,N);
    Legendre<double>::leg_der(D1,lval,gaussAbscissa,N,N);
    Legendre<double>::stiffLeg(stiff,gaussWeights,D1,N);

    for(outerLupe=0;outerLupe<N;++outerLupe)
        std::cout << gaussAbscissa[outerLupe] << ",";
    std::cout << gaussAbscissa[N] << std::endl;

    // Build an approximation to an ODE
    for(outerLupe=0;outerLupe<=N;++outerLupe)
        for(innerLupe=0;innerLupe<=N;++innerLupe)
        {
            jacobian[outerLupe][innerLupe] = stiff[outerLupe][innerLupe];
        }
    for(outerLupe=0;outerLupe<=N;++outerLupe)
    {
        jacobian[outerLupe][outerLupe] -= gaussWeights[outerLupe];
        baseFunc[outerLupe] = 1.0*gaussWeights[outerLupe];
    }

    if(!LU_Decomposition<double>::lu_decomp(jacobian,order,N+1))
    {
        std::cout << "The matrix is singular" << std::endl;
    }
    else
    {
        LU_Decomposition<double>::solve_lu(jacobian,deltaX,baseFunc,order,N+1);
        for(outerLupe=0;outerLupe<N;++outerLupe)
            std::cout << deltaX[outerLupe] << ",";
         std::cout << deltaX[N] << std::endl;
    }

    // Delete the arrays that have been allocated.
    ArrayUtils<double>::deltwotensor(lhsArray);
    ArrayUtils<double>::deltwotensor(secondDeriv);
    ArrayUtils<double>::deltwotensor(lval);
    ArrayUtils<double>::deltwotensor(stiff);
    ArrayUtils<double>::deltwotensor(D1);

    // Delete the vectors that have been allocated.
    ArrayUtils<double>::delonetensor(gaussWeights);
    ArrayUtils<double>::delonetensor(gaussAbscissa);

    // Delete the arrays and vectors used for the linear systems.
    ArrayUtils<double>::deltwotensor(jacobian);
    ArrayUtils<double>::delonetensor(baseFunc);
    ArrayUtils<double>::delonetensor(deltaX);
    ArrayUtils<int>::delonetensor(order);

    std::cout << "Done" << std::endl;
    return(0);
}
