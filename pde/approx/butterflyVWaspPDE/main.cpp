#include <iostream>

#include "util.h"
#include "legendre.h"


int main()
{
    std::cout << "Pre-processing starting" << std::endl;

    ArrayUtils<double> arrayCreator;
    Legendre<double> legendre;

    // Define the number of grid points to use.
    // Define the matrices that are used to create the system of equations.
    int N = 40;
    double **lhsArray    = arrayCreator.twotensor(N+1,N+1);
    double **secondDeriv = arrayCreator.twotensor(N+1,N+1);

    // Define the vectors that are used.
    // Includes the Gauss quadrature (abscissa and weights)
    double *gaussWeights  = arrayCreator.onetensor(N+1);
    double *gaussAbscissa = arrayCreator.onetensor(N+1);

    // Miscellaneous helper variables, for loops for example.
    int innerLupe;
    int outerLupe;

    // Define the Gauss quadrature.
    legendre.leg_quad(gaussAbscissa,gaussWeights,N);
    std::cout << "Abscissa: " << std::endl;
    for (innerLupe=0;innerLupe<=N;++innerLupe)
        std::cout << innerLupe << ": " << gaussAbscissa[innerLupe] << std::endl;

    std::cout << "Weights: " << std::endl;
    for (outerLupe=0;outerLupe<=N;++outerLupe)
        std::cout << outerLupe << ": " << gaussWeights[outerLupe] << std::endl;


    // Delete the arrays that have been allocated.
    arrayCreator.deltwotensor(lhsArray);
    arrayCreator.deltwotensor(secondDeriv);

    // Delete the vectors that have been allocated.
    arrayCreator.delonetensor(gaussWeights);
    arrayCreator.delonetensor(gaussAbscissa);
    return(0);
}
