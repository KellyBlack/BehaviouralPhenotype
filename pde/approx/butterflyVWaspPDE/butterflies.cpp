#include <iostream>

#include "butterflies.h"
#include "util.h"

Butterflies::Butterflies(int number, int sizeState) :
    PDESolver::PDESolver(number,sizeState)
{
    std::cout << "Butterflies begin" << std::endl;

    butterflies =nullptr;
    prevTimeStep =nullptr;

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
        integralSum   += gaussWeights[outerLupe]*0.5*dt*getG()*parameterDistribution(theta)*butterflies[outerLupe]/(c+butterflies[outerLupe]*parameterDistribution(theta));
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

//#define INITIAL_METHOD_ONE
void Butterflies::initializeButterflies()
{
    int number = getNumber();
    if(number>0)
    {
        if(butterflies==nullptr)
             butterflies  = ArrayUtils<double>::onetensor(number+1);
        if(prevTimeStep==nullptr)
             prevTimeStep = ArrayUtils<double>::onetensor(number+1);

        double integral1 = 0.0;
        double integral2 = 0.0;
        for(int lupe=0;lupe<=number;++lupe)
         {
             double theta = 0.5*(gaussAbscissa[lupe]+1.0);
             double steady = d*c/(parameterDistribution(theta)*(g-d));
             if(steady>1.0)
             {
                 butterflies[lupe]  = 1.0;
             }
             else if (steady>0.0)
             {
                 butterflies[lupe]  = steady;
             }
             else
             {
                 butterflies[lupe] = 0.01;
             }
             //integral1 += gaussWeights[lupe]*(1.0-butterflies[lupe])*(c+butterflies[lupe]*parameterDistribution(theta));
#ifdef INITIAL_METHOD_ONE
             integral1 += gaussWeights[lupe]*(c+butterflies[lupe]*parameterDistribution(theta));
             integral2 += gaussWeights[lupe]*(1.0-butterflies[lupe]);
#else
             integral1 += gaussWeights[lupe]*parameterDistribution(theta)*butterflies[lupe]/(c+butterflies[lupe]*parameterDistribution(theta));
             integral2 += gaussWeights[lupe]*parameterDistribution(theta)*butterflies[lupe]*(1.0-butterflies[lupe]);
#endif
         }
#ifdef INITIAL_METHOD_ONE
         butterflies[number+1] = integral2*integral1;
#else
         butterflies[number+1] = integral2/integral1;
#endif

         //double steady = d*c/(parameterDistribution(1.0)*(g-d));
         //butterflies[number+1] = (1.0-steady)*(c+steady*parameterDistribution(1.0));
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

void Butterflies::writeBinaryCurrentApprox(double &time,std::fstream &resultsFile)
{
    resultsFile.write(reinterpret_cast<char*>(&time),sizeof(double));
    resultsFile.write(reinterpret_cast<char*>(butterflies),static_cast<long>(N+2)*static_cast<long>(sizeof(double)));
}

double Butterflies::parameterDistribution(double theta)
{
    return(theta+a);
}
