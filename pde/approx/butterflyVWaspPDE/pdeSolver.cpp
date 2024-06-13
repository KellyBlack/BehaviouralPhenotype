#include <iostream>

#include "pdeSolver.h"

// Constructor for the PDE Solver class
PDESolver::PDESolver(int number, int numberState)
{
    // Set the requisite number of variables to keep track
    // of in the state of the PDE vector and the whole state
    // space as well.
    setNumber(number);
    setStateSize(numberState);

    // If the size of the state space is known allocate the
    // necessary space.
    if((number>0)&&(numberState>0))
        createArrays();

}

// Destructor for the PDE Solver class.
PDESolver::~PDESolver()
{
    // Just delete the vectors and arrays that have been allocated.
    deleteArrays();
}

// Mutator for the number variable.
void PDESolver::setNumber(int number)
{
    N = number;
}

// Accessor for the number variable
int PDESolver::getNumber()
{
    return(N);
}

// Mutator for the state space size.
void PDESolver::setStateSize(int number)
{
    stateSize = number;
}

// Accessor for the state space size.
int PDESolver::getStateSize()
{
    return(stateSize);
}

// Mutator for the time step variable.
void PDESolver::setDT(double value)
{
    dt = value;
}

// Accessor for the time step variable.
double PDESolver::getDT()
{
    return(dt);
}

// Routine to allocate the space used for the vectors and
// arrays associated with the class.
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

// Method to delete all the vectors and arrays used by the class.
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

// Method to initialize the values associated with the
// numerical scheme. In particular the quadrature values
// and the related stiffness and derivative matrices.
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


// Method to perform an LU decomposition on the current linear system,
// and then solve the decomposed system.
bool PDESolver::solveLinearizedSystem()
{
    LU_Decomposition<double> solver;
    if(solver.lu_decomp(jacobian,order,stateSize))
    {
        solver.solve_lu(jacobian,deltaX,baseFunc,order,stateSize);
        return(true);
    }
    return(false);
}

// Method to approximate the l2 norm of the update vector.
double PDESolver::normDelta()
{
    double value(0.0);
    int lupe;
    double *p = deltaX;
    for(lupe=0;lupe<stateSize;++lupe,++p)
    {
        value += (*p)*(*p);
    }
    return(value);
}

// Method to perform a single time step of the numerical scheme.
// Performs a Newton method to invert the nonlinear system.
int PDESolver::singleTimeStep(double maxNewtonDiffNorm,int maxNewtonSteps,bool printInfo)
{
    // Build the system and solve.
    calculateRHSTimeStepping();
    //copyCurrentStateToTemp();
    bool canInvert(true);
    double stepDeltaNorm = 0.0;
    int totalStepsPossible = maxNewtonSteps;
    do
    {
        // Perform the Newton steps to approximate the nonlinear
        // equations associated with the implicit system.
        buildJacobianTimeStepping();
        canInvert = solveLinearizedSystem();
        if(canInvert)
        {
            // Life is good. Perform the Newton step and update the current state vector.
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

// Method to approximate the steady state to the numerical approximation.
// Performs a Newton method to invert the nonlinear system.
int PDESolver::steadyStateApprox(double maxNewtonDiffNorm,int maxNewtonSteps,bool printInfo)
{
    // Build the system and solve.
    calculateRHSTimeStepping();
    bool canInvert(true);
    double stepDeltaNorm = 0.0;
    int totalStepsPossible = maxNewtonSteps;
    do
    {
        // Perform the Newton steps to approximate the nonlinear
        // equations associated with the implicit system.
        buildjacobianSteadyState();
        canInvert = solveLinearizedSystem();
        if(canInvert)
        {
            // Life is good. Perform the Newton step and update the current state vector.
            updateNewtonStep();
            if(printInfo)
                std::cout << "  step: " << totalStepsPossible - maxNewtonSteps << std::endl;
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

// Method to write the abscissa and weights associated with the
// Gauss-Labotto quadrature to a given file. The file is in a
// csv format and should be a text file.
void PDESolver::writeAbscissa(std::ofstream &resultsFile)
{
    for(int outerLupe=0;outerLupe<N;++outerLupe)
        resultsFile << gaussAbscissa[outerLupe] << ",";
    resultsFile << gaussAbscissa[N] << std::endl;

    for(int outerLupe=0;outerLupe<N;++outerLupe)
        resultsFile << gaussWeights[outerLupe] << ",";
    resultsFile << gaussWeights[N] << std::endl;

}

// Method to write the polynomial degree, the size of the state, and the
// abscissa and weights associated with the Gauss-Labotto quadrature to
// a binary file.
void PDESolver::writeBinaryHeader(std::fstream &resultsFile)
{
    resultsFile.write(reinterpret_cast<char*>(&N),sizeof(int));
    resultsFile.write(reinterpret_cast<char*>(&stateSize),sizeof(int));
    resultsFile.write(reinterpret_cast<char*>(gaussAbscissa),static_cast<long>(N+1)*static_cast<long>(sizeof(double)));
    resultsFile.write(reinterpret_cast<char*>(gaussWeights),static_cast<long>(N+1)*static_cast<long>(sizeof(double)));
}

