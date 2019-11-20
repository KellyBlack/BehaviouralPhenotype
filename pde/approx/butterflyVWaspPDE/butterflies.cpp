#include <iostream>

#include "butterflies.h"
#include "util.h"

Butterflies::Butterflies(int number, int sizeState) :
    PDESolver::PDESolver(number,sizeState)
{
    std::cout << "Butterflies begin" << std::endl;

    if(number>0)
        initializeButterflies();
}

Butterflies::~Butterflies()
{
    deleteButterflies();
    std::cout << "Butterflies End" << std::endl;
}


void Butterflies::buildJacobian()
{
    int outerLupe;
    int innerLupe;

    // Build an approximation to an ODE
    for(outerLupe=0;outerLupe<=N;++outerLupe)
    {
        baseFunc[outerLupe] = 0.0;
        for(innerLupe=0;innerLupe<=N;++innerLupe)
        {
            jacobian[outerLupe][innerLupe] = -0.5*dt*mu*stiff[outerLupe][innerLupe];
            baseFunc[outerLupe] -= stiff[outerLupe][innerLupe]*butterflies[innerLupe];
        }
        baseFunc[outerLupe] *= 0.5*dt*mu;
    }

    baseFunc[N+1] = butterflies[N+1]*(1.0+0.5*dt*getD())-rhs[N+1];
    double theta;
    double integralSum = 0.0;
    for(outerLupe=0;outerLupe<=N;++outerLupe)
    {
        theta = 0.5*(gaussAbscissa[outerLupe]+1.0);
        jacobian[outerLupe][outerLupe] += gaussWeights[outerLupe]*
                (
                    1.0-0.5*dt*(1.0-2.0*butterflies[outerLupe])*parameterDistribution(theta)
                    + 0.5*dt*butterflies[N+1]*parameterDistribution(theta)*c/((c+butterflies[outerLupe]*parameterDistribution(theta))*(c+butterflies[outerLupe]*parameterDistribution(theta)))
                 );
        baseFunc[outerLupe] += gaussWeights[outerLupe]*(
                    butterflies[outerLupe] -
                    0.5*dt*parameterDistribution(theta)*butterflies[outerLupe]*(1.0-butterflies[outerLupe]) +
                    0.5*dt*butterflies[N+1]*parameterDistribution(theta)*butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta)))
                - rhs[outerLupe];

        jacobian[N+1][outerLupe] =
                gaussWeights[outerLupe]*0.5*dt*getG()*butterflies[N+1]*parameterDistribution(theta)*c/((c+butterflies[outerLupe]*parameterDistribution(theta))*(c+butterflies[outerLupe]*parameterDistribution(theta)));
        baseFunc[N+1] -= gaussWeights[outerLupe]*0.5*dt*getG()*butterflies[N+1]*parameterDistribution(theta)*butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta));
        integralSum += gaussWeights[outerLupe]*0.5*dt*getG()*parameterDistribution(theta)*butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta));
    }
    //theta = 0.5*(gaussAbscissa[N]+1.0);
    jacobian[N+1][N+1] = 1.0+0.5*dt*getD() - integralSum;

}

void Butterflies::updateNewtonStep()
{
    int lupe;
    for(lupe=0;lupe<=N+1;++lupe)
    {
        butterflies[lupe] -= deltaX[lupe];
        //std::cout << lupe << ": " << deltaX[lupe] << "/" << butterflies[lupe] << std::endl;
    }
}

void Butterflies::calculateRHS()
{
    int outerLupe;
    int innerLupe;

    rhs[N+1] = butterflies[N+1]*(1.0-dt*0.5*getD());
    for(outerLupe=0;outerLupe<=N;++outerLupe)
    {
        double theta = 0.5*(gaussAbscissa[outerLupe]+1.0);
        rhs[outerLupe] = 0.0;
        for(innerLupe=0;innerLupe<=N;++innerLupe)
        {
            rhs[outerLupe] += stiff[outerLupe][innerLupe]*butterflies[innerLupe];
        }
        rhs[outerLupe] = 0.5*dt*mu*rhs[outerLupe] +
                gaussWeights[outerLupe]*(
                    butterflies[outerLupe] +
                    0.5*dt*parameterDistribution(theta)*butterflies[outerLupe]*(1.0-butterflies[outerLupe]) -
                    0.5*dt*parameterDistribution(theta)*butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta))
                    );
        rhs[N+1] += gaussWeights[outerLupe]*0.5*dt*getG()*butterflies[N+1]*parameterDistribution(theta)*butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta));
    }
}

void Butterflies::copyCurrentStateToTemp()
{
    int number = getNumber();
    for(int lupe=0;lupe<=number;++lupe)
        prevTimeStep[lupe] = butterflies[lupe];
}

void Butterflies::initializeButterflies()
{
    int number = getNumber();
    if(butterflies!=nullptr)
        ArrayUtils<double>::delonetensor(butterflies);

    if(prevTimeStep!=nullptr)
        ArrayUtils<double>::delonetensor(prevTimeStep);

    if(number>0)
    {
         butterflies  = ArrayUtils<double>::onetensor(number+1);
         prevTimeStep = ArrayUtils<double>::onetensor(number+1);
         for(int lupe=0;lupe<=number;++lupe)
         {
             butterflies[lupe]  = 0.9;
         }
         butterflies[number+1] = 0.5;
         copyCurrentStateToTemp();
    }
}

void Butterflies::deleteButterflies()
{
    if(butterflies!=nullptr)
        ArrayUtils<double>::delonetensor(butterflies);
}

void Butterflies::writeCurrentApprox(double time, std::ofstream &resultsFile)
{
    resultsFile << time;
    for(int outerLupe=0;outerLupe<=N;++outerLupe)
        resultsFile << "," << butterflies[outerLupe];
    resultsFile << std::endl;

}

double Butterflies::parameterDistribution(double theta)
{
    return((theta+1.0)*0.5);
}
