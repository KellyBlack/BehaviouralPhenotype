#include "butterflies.h"

butterflies::butterflies(int number)
{
    setNumber(number);
    if(number>0)
        createArrays();

}

butterflies::~butterflies()
{
    deleteArrays();
}

void butterflies::setNumber(int number)
{
    N = number;
}

int butterflies::getNumber()
{
    return(N);
}

void butterflies::createArrays()
{
    lhsArray    = ArrayUtils<double>::twotensor(N+1,N+1);
    secondDeriv = ArrayUtils<double>::twotensor(N+1,N+1);
    lval        = ArrayUtils<double>::twotensor(N+1,N+1);
    D1          = ArrayUtils<double>::twotensor(N+1,N+1);
    stiff       = ArrayUtils<double>::twotensor(N+1,N+1);

    // Define the vectors that are used.
    // Includes the Gauss quadrature (abscissa and weights)
    gaussWeights  = ArrayUtils<double>::onetensor(N+1);
    gaussAbscissa = ArrayUtils<double>::onetensor(N+1);

    // Define the matrices and vectors used to solve the linear systems.
    jacobian = ArrayUtils<double>::twotensor(N+1,N+1);
    baseFunc  = ArrayUtils<double>::onetensor(N+1);
    deltaX  = ArrayUtils<double>::onetensor(N+1);
    order     = ArrayUtils<int>::onetensor(N+1);
}

void butterflies::deleteArrays()
{
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

}
