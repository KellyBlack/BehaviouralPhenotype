#include <iostream>
#include <cmath>

#include "butterflies.h"
#include "util.h"

// Constructor for the Butterfly class.
Butterflies::Butterflies(int number, int sizeState) :
    PDESolver::PDESolver(number,sizeState)
{
    // If a size was passed, so go ahead and allocate the
    // matrices and set up the initial condition.
    if(getNumber()>0)
        initializeButterflies();
}

// Destructor for the class.
Butterflies::~Butterflies()
{
    // De-allocate the arrays and vectors that are in use.
    deleteButterflies();
}


// Routine to set up the entries in the Jacobian for the
// full nonlinear implicit time step.
void Butterflies::buildJacobian()
{
    // integers used in various loops.
    int outerLupe;
    int innerLupe;

    // Build an approximation to an ODE
    // Set up the pointers that are used to step through the
    // columns of the matrix. (Assumes the arrays are in
    // row major order.)
    double *j;
    double *s;
    double *b = baseFunc;

    // Add the coefficients of the stiffness matrix for the
    // whole Jacobian assocated with the butterfly state vector.
    for(outerLupe=0;outerLupe<=N;++outerLupe)
    {
        // Initialize the base function to zero.
        *b = 0.0;
        j = jacobian[outerLupe];
        s = stiff[outerLupe];

        // For each column in the current row set the entries in the
        // Jacobian to be the associated entries in the stiffness matrix.
        // Add up the dot product with the row in the stiffness matrix
        // and the butterfly state vector for the function base function.
        for(innerLupe=0;innerLupe<=N;++innerLupe)
        {
            *j++ = -0.5*dt*mu*(*s)*4.0;                // 2.0 = 0.5*4.0 (0.5 from time stepping 4.0 from spatial mapping)
            *b -= (*s++)*butterflies[innerLupe];
        }
        // Multiply the stiffness matrix product by the diffusion term and
        // the time step associated with the implicit scheme.
        *b++ *= 0.5*dt*mu*4.0;     // 2.0 = 0.5*4.0 (0.5 from time stepping 4.0 from spatial mapping)
    }

    // Now go through and add the terms associated with the diagonal entries
    // of the Jacobian as well as the terms associated with the equation
    // for the wasps.

    // Initialize the base function term associated with the wasps.
    baseFunc[N+1] = butterflies[N+1]*(1.0+0.5*dt*getD())-rhs[N+1];
    double theta;
    double integralSum = 0.0;
    for(outerLupe=0;outerLupe<=N;++outerLupe)
    {
        // Add all of the derivatives associated with the diagonal entries of the Jacobian.
        theta = 0.5*(gaussAbscissa[outerLupe]+1.0);
        jacobian[outerLupe][outerLupe] += gaussWeights[outerLupe]*
                (
                    1.0-0.5*dt*(1.0-2.0*butterflies[outerLupe])*parameterDistribution(theta)
                    + 0.5*dt*butterflies[N+1]*parameterDistribution(theta)*c/((c+butterflies[outerLupe]*parameterDistribution(theta))*(c+butterflies[outerLupe]*parameterDistribution(theta)))
                 );

        // Add all of the function values associated with the diagonal entries to the base function.
        baseFunc[outerLupe] += gaussWeights[outerLupe]*(
                    butterflies[outerLupe] -
                    0.5*dt*parameterDistribution(theta)*butterflies[outerLupe]*(1.0-butterflies[outerLupe]) +
                    0.5*dt*butterflies[N+1]*parameterDistribution(theta)*butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta)))
                - rhs[outerLupe];

        // Add all of the remaining nonlinear terms to the bottom row of the Jacobian.
        // These are the partial derivatives associated with the integral of the butterfly state space.
        // 0.5 from time stepping and 0.5 from mapping the integral to theta in [0,1] from xi in [-1,1].
        jacobian[N+1][outerLupe] =
                0.5*gaussWeights[outerLupe]*0.5*dt*getG()*butterflies[N+1]*parameterDistribution(theta)*c/((c+butterflies[outerLupe]*parameterDistribution(theta))*(c+butterflies[outerLupe]*parameterDistribution(theta)));

        // Add the integral terms to the base function for the wasps. Also keep track of the partial derivative of
        // the integral with respect to the wasps.
        // 0.5 from time stepping and 0.5 from mapping the integral to theta in [0,1] from xi in [-1,1].
        baseFunc[N+1] -= 0.5*gaussWeights[outerLupe]*0.5*dt*getG()*butterflies[N+1]*parameterDistribution(theta)*butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta));
        integralSum   += 0.5*gaussWeights[outerLupe]*0.5*dt*getG()*parameterDistribution(theta)*butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta));
    }

    // Finally update the Jacobian for the partial derivative of the wasps with respect to the wasps.
    jacobian[N+1][N+1] = 1.0+0.5*dt*getD() - integralSum;

}

// Method to add the deltaX from the Newton method solve to the current state space.
void Butterflies::updateNewtonStep()
{
    int lupe;
    double *b = butterflies;
    double *dx = deltaX;
    for(lupe=0;lupe<=N+1;++lupe)
    {
        *b++ -= *dx++;
    }
}

// Method to calculate the function value associated with the
// previous time step. These are the known values of the time
// stepping equation.
void Butterflies::calculateRHS()
{
    // loop variables.
    int outerLupe;
    int innerLupe;

    // Initialize the values associated with the linear terms in the equation for the wasps.
    rhs[N+1] = butterflies[N+1]*(1.0-dt*0.5*getD());

    // Initialize the pointer to the rhs vector.
    double *r = rhs;
    double *s;
    double *b;
    for(outerLupe=0;outerLupe<=N;++outerLupe)
    {
        // Initialize the value in the RHS vector to be zero.
        // Later we will add the terms from the function evaluation.
        double theta = 0.5*(gaussAbscissa[outerLupe]+1.0);
        *r = 0.0;

        // Initialize the pointers to the current row of the stiffness matrix
        // and the starting of the state vector for the butterflies.
        s = stiff[outerLupe];
        b = butterflies;
        for(innerLupe=0;innerLupe<=N;++innerLupe)
        {
            // Add the terms to the RHS for the dot product of the
            // current row in the stiffness matrix and the butterfly
            // state vector.
            // 4.0 comes from mapping theta in [0,1] from xi in [-1,1].
            *r += (*s++)*(*b++)*4.0;
        }

        // Add all of the terms associated with the remaining terms in the butterfly
        // equation. (All of the terms not associated with the stiffness matrix.)
        *r++ = 0.5*dt*mu*rhs[outerLupe] +
                gaussWeights[outerLupe]*(
                    butterflies[outerLupe] +
                    0.5*dt*parameterDistribution(theta)*butterflies[outerLupe]*(1.0-butterflies[outerLupe]) -
                    0.5*dt*parameterDistribution(theta)*butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta))
                    );

        // Add the terms associated with the integral in the wasp equation.
        // 0.5comes from mapping theta in [0,1] from xi in [-1,1].
        rhs[N+1] += 0.5*gaussWeights[outerLupe]*0.5*dt*getG()*butterflies[N+1]*parameterDistribution(theta)*
                butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta));
    }
}

// Method to copy the current state vector to the prevTimeStep vector.
// Is this really needed???? Does not seem to be necessary?
void Butterflies::copyCurrentStateToTemp()
{
    int number = getStateSize();
    double *p = prevTimeStep;
    double *b = butterflies;
    for(int lupe=0;lupe<number;++lupe)
        *p++ = *b++;
}

// Method to copy the current state vector to a vector that
// is passed in as the argument.
void Butterflies::copyCurrentState(double *ptr)
{
    double *dest = ptr;
    double *source = butterflies;
    for(int lupe=0;lupe<=getStateSize();++lupe)
        *dest++ = *source++;
}

// Method to approximate the total butterfly population by
// integrating the current state.
double Butterflies::totalButterflyPopulation()
{
    double integral = 0.0;
    double *b = butterflies;
    double *gw = gaussWeights;
    for(int lupe=0;lupe<=N;++lupe)
    {
        integral += (*b++)*(*gw++);
    }
    return(integral*0.5);
}

// Method to return the current wasp density
double Butterflies::waspPopulation()
{
    return(butterflies[getStateSize()-1]);
}

// Method to allocate space for the state vector and also set the initial
// conditions for the butterflies and the wasps.
void Butterflies::initializeButterflies()
{
    int number = getNumber();
    if(number>0)
    {
        // This object know how big it should be. Allocate the space
        // if it has not already been done.
        ArrayUtils <double> arrays;
        if(butterflies==nullptr)
             butterflies  = arrays.onetensor(getStateSize());
        if(prevTimeStep==nullptr)
             prevTimeStep = arrays.onetensor(getStateSize());

        // First set the initial profile for the butterflies. It will be normalized
        // later so its average value is something more reasonable.
        double butterflyIntegral = 0.0;
        for(int lupe=0;lupe<=number;++lupe)
         {
             // Use the theoretical steady state predicted by the ODE.
             double theta = 0.5*(gaussAbscissa[lupe]+1.0);
             double steady = d*c/(parameterDistribution(theta)*(g-d));
             /*
             if(steady>1.0)
             {
                 butterflies[lupe]  = 1.0;
             } else
             */
             if (steady>0.0)
             {
                 // The predicted steady state is positive so assume it is okay.
                 butterflies[lupe]  = steady;
             }
             else
             {
                 // The predicted steady state was negative so just make it
                 // something close to zero.
                 butterflies[lupe] = 0.01;
             }

             // Keep track of the integral of the butterflies to estimate
             // the total number of butterflies in the population.
             butterflyIntegral += gaussWeights[lupe]*butterflies[lupe];
         }

//#define INITIAL_METHOD_ONE
        // Normalize the butterfly density so that the total number
        // is a reasonable overall value. Reuse the butterflyIntegral
        // variable to be the new scaling factor.
        double integral1 = 0.0;
        double integral2 = 0.0;
        butterflyIntegral = d*c/(parameterDistribution(0.0)*(g-d))/butterflyIntegral;
        // Now normalize the butterflies and integrate the butterfly DE assuming
        // steady state to get the initial approximation for the wasps.
        for(int lupe=0;lupe<=number;++lupe)
         {
            // Go through and normalize every value in the state vector associated
            // with the butterflies.
            double theta = 0.5*(gaussAbscissa[lupe]+1.0);
            butterflies[lupe] *= butterflyIntegral;
#ifdef INITIAL_METHOD_ONE
            // Integrate the terms in the steady state equation for each butterfly value.
             integral1 += gaussWeights[lupe]*(c+butterflies[lupe]*parameterDistribution(theta));
             integral2 += gaussWeights[lupe]*(1.0-butterflies[lupe]);
#else
            // Integrate the terms in the full butterfly, time varying, equation assuming that the time derivative is zero.
             integral1 += gaussWeights[lupe]*parameterDistribution(theta)*butterflies[lupe]/(c+butterflies[lupe]*parameterDistribution(theta));
             integral2 += gaussWeights[lupe]*parameterDistribution(theta)*butterflies[lupe]*fabs((1.0-butterflies[lupe]));
#endif
         }
#ifdef INITIAL_METHOD_ONE
        // Make an estimate for the number of wasps based on the steady state equation.
         butterflies[number+1] = integral2*integral1;
#else
        // Make an estimate for the number of wasps based on the full PDE assuming the time derivative is zero.
         butterflies[number+1] = integral2/integral1;
#endif

         copyCurrentStateToTemp();
    }
}

// Method to allocate space for the state vector and also set the initial
// conditions for the butterflies and the wasps given an initial value.
void Butterflies::initializeButterflies(double *initial)
{
    int number = getNumber();
    if(number>0)
    {
        // This object know how big it should be. Allocate the space
        // if it has not already been done.
        ArrayUtils <double> arrays;
        if(butterflies==nullptr)
             butterflies  = arrays.onetensor(getStateSize());
        if(prevTimeStep==nullptr)
             prevTimeStep = arrays.onetensor(getStateSize());

        // First set the initial profile for the butterflies. It will be normalized
        // later so its average value is something more reasonable.
        for(int lupe=0;lupe<=number;++lupe)
         {
            butterflies[lupe] = initial[lupe];
         }
         butterflies[number+1] = initial[number+1];
         copyCurrentStateToTemp();
    }
}


// Method to set the initial condition to a Gaussian profile with a given
// center and variance.
void Butterflies::initializeButterfliesGaussian(double center,double variance)
{
    int number = getNumber();
    if(number>0)
    {
        // Allocate the space if the size of the state space is known.
        ArrayUtils <double> arrays;
        if(butterflies==nullptr)
             butterflies  = arrays.onetensor(getStateSize());
        if(prevTimeStep==nullptr)
             prevTimeStep = arrays.onetensor(getStateSize());

        // First set the initial profile for the butterflies. It is a
        // simple Gaussian with the given center and variance.
        double butterflyIntegral = 0.0;
        for(int lupe=0;lupe<=number;++lupe)
         {
            double theta = 0.5*(gaussAbscissa[lupe]+1.0);
            butterflies[lupe]  = exp(-(theta-center)*(theta-center)/variance);
            butterflyIntegral += gaussWeights[lupe]*butterflies[lupe];
         }

        double steady = d*c/(parameterDistribution(0.5)*(g-d));
        butterflies[number+1] = fabs(1.0-steady)*(c+steady*parameterDistribution(0.0));
        steady /= (0.5*butterflyIntegral);
        /*
        for(int lupe=0;lupe<=number;++lupe)
         {
            butterflies[lupe] *= steady;
         }
         */
         copyCurrentStateToTemp();
    }
}

// Method to set the initial condition to a constant with the
// steady state at the given value of theta.
void Butterflies::initializeButterfliesConstant(double theta)
{
    int number = getNumber();
    if(number>0)
    {
        // Allocate the space if the size of the state space is known.
        ArrayUtils <double> arrays;
        if(butterflies==nullptr)
             butterflies  = arrays.onetensor(getStateSize());
        if(prevTimeStep==nullptr)
             prevTimeStep = arrays.onetensor(getStateSize());

        // First set the initial profile for the butterflies. It is a
        // constant at the given value of theta.
        double steady = d*c/(parameterDistribution(theta)*(g-d));
        for(int lupe=0;lupe<=number;++lupe)
         {
            butterflies[lupe]  = steady;
         }
        butterflies[number+1] = fabs(1.0-steady)*(c+steady*parameterDistribution(theta));
         copyCurrentStateToTemp();
    }
}

// Method to write the relevant parameter values to
// a given binary file.
void Butterflies::writeParameters(std::fstream &resultsFile)
{
    resultsFile.write(reinterpret_cast<char*>(&mu),sizeof(double));
    resultsFile.write(reinterpret_cast<char*>(&c),sizeof(double));
    resultsFile.write(reinterpret_cast<char*>(&g),sizeof(double));
    resultsFile.write(reinterpret_cast<char*>(&d),sizeof(double));
    resultsFile.write(reinterpret_cast<char*>(&m),sizeof(double));
}

// Method to delete the state vector and the
// vector used to keep track of the previous time step.
void Butterflies::deleteButterflies()
{
    ArrayUtils<double> arrays;
    if(butterflies!=nullptr)
        arrays.delonetensor(butterflies);
    butterflies = nullptr;

    if(prevTimeStep!=nullptr)
        arrays.delonetensor(prevTimeStep);
    prevTimeStep = nullptr;

}

// Method to write the current state vector to a csv file.
// The file should be opened as a standard text file.
void Butterflies::writeCurrentApprox(double time, std::ofstream &resultsFile)
{
    resultsFile << time;
    for(int outerLupe=0;outerLupe<=N;++outerLupe)
        resultsFile << "," << butterflies[outerLupe];
    resultsFile << "," << butterflies[N+1] << std::endl; // write out the wasp density!

}

// Method to write the current state vector to a binary file.
void Butterflies::writeBinaryCurrentApprox(double &time,std::fstream &resultsFile)
{
    // First approximate the total butterfly population.
    double integral = totalButterflyPopulation();

    // Write the time, the total population, and then the full state vector.
    resultsFile.write(reinterpret_cast<char*>(&time),sizeof(double));
    resultsFile.write(reinterpret_cast<char*>(&integral),sizeof(double));
    resultsFile.write(reinterpret_cast<char*>(butterflies),static_cast<long>(N+2)*static_cast<long>(sizeof(double)));
}

// Method to calculate the value of p(theta).
double Butterflies::parameterDistribution(double theta)
{
    return(m*theta+1.0);
}
