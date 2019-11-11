#include <iostream>

#include "pdeSolver.h"

PDESolver::PDESolver(int number)
{
    std::cout << "PDE Solver start" << std::endl;
    setNumber(number);
    if(number>0)
        createArrays();

}

PDESolver::~PDESolver()
{
    std::cout << "PDE Solver end" << std::endl;
    deleteArrays();
}

void PDESolver::setNumber(int number)
{
    N = number;
}

int PDESolver::getNumber()
{
    return(N);
}

void PDESolver::createArrays()
{
    // Define the matrices used to store necessary values for
    // the Legendre polynomials.
    lval        = ArrayUtils<double>::twotensor(N+1,N+1);
    D1          = ArrayUtils<double>::twotensor(N+1,N+1);
    stiff       = ArrayUtils<double>::twotensor(N+1,N+1);

    // Define the vectors that are used.
    // Includes the Gauss quadrature (abscissa and weights)
    gaussWeights  = ArrayUtils<double>::onetensor(N+1);
    gaussAbscissa = ArrayUtils<double>::onetensor(N+1);

    // Define the matrices and vectors used to solve the linear systems.
    jacobian = ArrayUtils<double>::twotensor(N+1,N+1);
    baseFunc = ArrayUtils<double>::onetensor(N+1);
    deltaX   = ArrayUtils<double>::onetensor(N+1);
    order    = ArrayUtils<int>::onetensor(N+1);
}

void PDESolver::deleteArrays()
{
    // Delete the arrays that have been allocated.
    if(lval!=nullptr)
        ArrayUtils<double>::deltwotensor(lval);
    if(stiff!=nullptr)
        ArrayUtils<double>::deltwotensor(stiff);
    if(D1!=nullptr)
        ArrayUtils<double>::deltwotensor(D1);

    // Delete the vectors that have been allocated.
    if(gaussWeights!=nullptr)
        ArrayUtils<double>::delonetensor(gaussWeights);
    if(gaussAbscissa!=nullptr)
        ArrayUtils<double>::delonetensor(gaussAbscissa);

    // Delete the arrays and vectors used for the linear systems.
    if(jacobian!=nullptr)
        ArrayUtils<double>::deltwotensor(jacobian);
    if(baseFunc!=nullptr)
        ArrayUtils<double>::delonetensor(baseFunc);
    if(deltaX!=nullptr)
        ArrayUtils<double>::delonetensor(deltaX);
    if(order!=nullptr)
        ArrayUtils<int>::delonetensor(order);

}

void PDESolver::initializeLegendreParams()
{
    // Define the Gauss quadrature.
    Legendre<double>::leg_quad(gaussAbscissa,gaussWeights,N);

    // Define the stiffness and mass matrices for the Legendre collocation method.
    Legendre<double>::leg_val(lval,gaussAbscissa,N,N);
    Legendre<double>::leg_der(D1,lval,gaussAbscissa,N,N);
    Legendre<double>::stiffLeg(stiff,gaussWeights,D1,N);

}


bool PDESolver::solveLinearizedSystem()
{
    if(LU_Decomposition<double>::lu_decomp(jacobian,order,N+1))
    {
        LU_Decomposition<double>::solve_lu(jacobian,deltaX,baseFunc,order,N+1);
        return(true);
    }
    return(false);
}

void PDESolver::writeAbscissa(std::ofstream &resultsFile)
{
    for(int outerLupe=0;outerLupe<N;++outerLupe)
        resultsFile << gaussAbscissa[outerLupe] << ",";
    resultsFile << gaussAbscissa[N] << std::endl;
}

void PDESolver::writeCurrentApprox(std::ofstream &resultsFile)
{
    for(int outerLupe=0;outerLupe<N;++outerLupe)
        resultsFile << deltaX[outerLupe] << ",";
     resultsFile << deltaX[N] << std::endl;
}
