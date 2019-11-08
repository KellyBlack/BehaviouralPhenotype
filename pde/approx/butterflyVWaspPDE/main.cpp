#include <iostream>
#include <cmath>

#include "util.h"
#include "legendre.h"


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

    // Miscellaneous helper variables, for loops for example.
    int innerLupe;
    int outerLupe;

    // Define the Gauss quadrature.
    Legendre<double>::leg_quad(gaussAbscissa,gaussWeights,N);
    std::cout << "Abscissa: " << std::endl;
    for (innerLupe=0;innerLupe<=N;++innerLupe)
        std::cout << innerLupe << ": " << gaussAbscissa[innerLupe] << std::endl;

    std::cout << "Weights: " << std::endl;
    for (outerLupe=0;outerLupe<=N;++outerLupe)
        std::cout << outerLupe << ": " << gaussWeights[outerLupe] << std::endl;

    // Define the stiffness and mass matrices for the Legendre collocation method.
    std::cout << "Pre-processing" << std::endl;
    Legendre<double>::leg_val(lval,gaussAbscissa,N,N);
    Legendre<double>::leg_der(D1,lval,gaussAbscissa,N,N);
    Legendre<double>::stiffLeg(stiff,gaussWeights,D1,N);

    // Test the D1 matrix
    std::cout << "Testing" << std::endl;
    //double x2 = 0.0;
    double sum;
#define POWER 20
    for(outerLupe=0;outerLupe<=N;++outerLupe)
    {
        sum = 0.0;
        for(innerLupe=0;innerLupe<=N;++innerLupe)
        {
            sum += D1[outerLupe][innerLupe]*pow(gaussAbscissa[innerLupe],static_cast<double>(POWER));
                   //gaussWeights[innerLupe]*pow(gaussAbscissa[innerLupe],static_cast<double>(outerLupe));
        }
        std::cout << "Deriv: " << sum << " - " << static_cast<double>(POWER)*pow(gaussAbscissa[outerLupe],static_cast<double>(POWER-1)) << std::endl;
        //x2 += 1.0;
        //std::cout << "integral: " << sum << " - " << (pow(1.0,static_cast<double>(outerLupe+1))-pow(-1.0,static_cast<double>(outerLupe+1)))/x2 << std::endl;
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

    std::cout << "Done" << std::endl;
    return(0);
}
