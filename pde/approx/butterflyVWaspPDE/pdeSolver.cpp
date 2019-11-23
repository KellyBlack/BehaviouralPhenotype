#include <iostream>

#include "pdeSolver.h"

PDESolver::PDESolver(int number, int numberState)
{
    std::cout << "PDE Solver start" << std::endl;
    setNumber(number);
    setStateSize(numberState);
    if((number>0)&&(numberState>0))
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

void PDESolver::setStateSize(int number)
{
    stateSize = number;
}

int PDESolver::getStateSize()
{
    return(stateSize);
}


void PDESolver::setDT(double value)
{
    dt = value;
}

double PDESolver::getDT()
{
    return(dt);
}

void PDESolver::createArrays()
{
    // First delete all of the current arrays.
    deleteArrays();

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
    jacobian = ArrayUtils<double>::twotensor(stateSize,stateSize);
    baseFunc = ArrayUtils<double>::onetensor(stateSize);
    rhs      = ArrayUtils<double>::onetensor(stateSize);
    deltaX   = ArrayUtils<double>::onetensor(stateSize);
    order    = ArrayUtils<int>::onetensor(stateSize);
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
    if(rhs!=nullptr)
        ArrayUtils<double>::delonetensor(rhs);
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
    if(LU_Decomposition<double>::lu_decomp(jacobian,order,stateSize))
    {
        LU_Decomposition<double>::solve_lu(jacobian,deltaX,baseFunc,order,stateSize);
        return(true);
    }
    return(false);
}

double PDESolver::normDelta()
{
    double value(0.0);
    int lupe;
    for(lupe=0;lupe<stateSize;++lupe)
        value += deltaX[lupe]*deltaX[lupe];
    return(value);
}

void PDESolver::writeAbscissa(std::ofstream &resultsFile)
{
    for(int outerLupe=0;outerLupe<N;++outerLupe)
        resultsFile << gaussAbscissa[outerLupe] << ",";
    resultsFile << gaussAbscissa[N] << std::endl;

    for(int outerLupe=0;outerLupe<N;++outerLupe)
        resultsFile << gaussWeights[outerLupe] << ",";
    resultsFile << gaussWeights[N] << std::endl;

}

void PDESolver::writeBinaryHeader(std::fstream &resultsFile)
{
    resultsFile.write(reinterpret_cast<char*>(&N),sizeof(int));
    resultsFile.write(reinterpret_cast<char*>(&stateSize),sizeof(int));
    resultsFile.write(reinterpret_cast<char*>(gaussAbscissa),static_cast<long>(N+1)*static_cast<long>(sizeof(double)));
    resultsFile.write(reinterpret_cast<char*>(gaussWeights),static_cast<long>(N+1)*static_cast<long>(sizeof(double)));
}

