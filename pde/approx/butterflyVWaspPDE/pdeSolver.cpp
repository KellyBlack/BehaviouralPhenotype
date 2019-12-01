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
    ArrayUtils<double> arrays;
    ArrayUtils<int> arraysInt;
    // First delete all of the current arrays.
    deleteArrays();

    // Define the matrices used to store necessary values for
    // the Legendre polynomials.
    lval        = arrays.twotensor(N+1,N+1);
    D1          = arrays.twotensor(N+1,N+1);
    stiff       = arrays.twotensor(N+1,N+1);

    // Define the vectors that are used.
    // Includes the Gauss quadrature (abscissa and weights)
    gaussWeights  = arrays.onetensor(N+1);
    gaussAbscissa = arrays.onetensor(N+1);

    // Define the matrices and vectors used to solve the linear systems.
    jacobian = arrays.twotensor(stateSize,stateSize);
    baseFunc = arrays.onetensor(stateSize);
    rhs      = arrays.onetensor(stateSize);
    deltaX   = arrays.onetensor(stateSize);
    order    = arraysInt.onetensor(stateSize);
}

void PDESolver::deleteArrays()
{
    ArrayUtils<double> arrays;
    ArrayUtils<int> arraysInt;

    // Delete the arrays that have been allocated.
    if(lval!=nullptr)
    {
        arrays.deltwotensor(lval);
        lval = nullptr;
    }

    if(stiff!=nullptr)
    {
        arrays.deltwotensor(stiff);
        stiff  = nullptr;
    }

    if(D1!=nullptr)
    {
        arrays.deltwotensor(D1);
        D1  = nullptr;
    }

    // Delete the vectors that have been allocated.
    if(gaussWeights!=nullptr)
    {
        arrays.delonetensor(gaussWeights);
        gaussWeights = nullptr;
    }

    if(gaussAbscissa!=nullptr)
    {
        arrays.delonetensor(gaussAbscissa);
        gaussAbscissa = nullptr;
    }

    // Delete the arrays and vectors used for the linear systems.
    if(jacobian!=nullptr)
    {
        arrays.deltwotensor(jacobian);
        jacobian  = nullptr;
    }

    if(baseFunc!=nullptr)
    {
        arrays.delonetensor(baseFunc);
        baseFunc = nullptr;
    }

    if(rhs!=nullptr)
    {
        arrays.delonetensor(rhs);
        rhs = nullptr;
    }

    if(deltaX!=nullptr)
    {
        arrays.delonetensor(deltaX);
        deltaX = nullptr;
    }

    if(order!=nullptr)
    {
        arraysInt.delonetensor(order);
        order = nullptr;
    }
}

void PDESolver::initializeLegendreParams()
{
    if(getNumber()>0)
    {
        Legendre<double> legendre;
        // Define the Gauss quadrature.
        legendre.leg_quad(gaussAbscissa,gaussWeights,N);

        // Define the stiffness and mass matrices for the Legendre collocation method.
        legendre.leg_val(lval,gaussAbscissa,N,N);
        legendre.leg_der(D1,lval,gaussAbscissa,N,N);
        legendre.stiffLeg(stiff,gaussWeights,D1,N);
    }

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

int PDESolver::singleTimeStep(double maxNewtonDiffNorm,int maxNewtonSteps,bool printInfo)
{
    // Build the system and solve.
    calculateRHS();
    copyCurrentStateToTemp();
    bool canInvert(true);
    double stepDeltaNorm = 0.0;
    int totalStepsPossible = maxNewtonSteps;
    do
    {
        // Perform the Newton steps to approximate the nonlinear
        // equations associated with the implicit system.
        buildJacobian();
        canInvert = solveLinearizedSystem();
        if(canInvert)
        {
            updateNewtonStep();
            if(printInfo)
                std::cout << "  step: " << totalStepsPossible - maxNewtonSteps << ", ";
            stepDeltaNorm = normDelta();
        }
        else
        {
            std::cout << "  System not invertible" << std::endl;
            return(-1);
        }

    } while((stepDeltaNorm>maxNewtonDiffNorm) && canInvert && (maxNewtonSteps-- > 0));

    if(printInfo)
        std::cout << std::endl;

    return(totalStepsPossible - maxNewtonSteps);

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

