#include <iostream>

#include "butterflies.h"
#include "util.h"

Butterflies::Butterflies(int number) : PDESolver::PDESolver(number)
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
            jacobian[outerLupe][innerLupe] = mu*stiff[outerLupe][innerLupe];
            baseFunc[outerLupe] -= mu*stiff[outerLupe][innerLupe]*butterflies[innerLupe];
        }
    }

    for(outerLupe=0;outerLupe<=N;++outerLupe)
    {
        jacobian[outerLupe][outerLupe] += gaussWeights[outerLupe]*(1.0+(1.0-2.0*butterflies[outerLupe]));
        baseFunc[outerLupe] -= gaussWeights[outerLupe]*(
                    butterflies[outerLupe]*(1.0-butterflies[outerLupe]) +
                    prevTimeStep[outerLupe] +
                    butterflies[outerLupe]);
    }

}

void Butterflies::updateNewtonStep()
{
    int lupe;
    for(lupe=0;lupe<=N;++lupe)
        butterflies[lupe] += deltaX[lupe];
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
             butterflies[lupe]  = 1.0;
             prevTimeStep[lupe] = butterflies[lupe];
         }
    }
}

void Butterflies::deleteButterflies()
{
    if(butterflies!=nullptr)
        ArrayUtils<double>::delonetensor(butterflies);
}
