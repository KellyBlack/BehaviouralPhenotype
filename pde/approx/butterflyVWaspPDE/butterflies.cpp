#include <iostream>

#include "butterflies.h"

Butterflies::Butterflies(int number) : PDESolver::PDESolver(number)
{
    std::cout << "Butterflies begin" << std::endl;
}

Butterflies::~Butterflies()
{
    std::cout << "Butterflies End" << std::endl;
}


void Butterflies::buildJacobian()
{
    int outerLupe;
    int innerLupe;

    // Build an approximation to an ODE
    for(outerLupe=0;outerLupe<=N;++outerLupe)
        for(innerLupe=0;innerLupe<=N;++innerLupe)
        {
            jacobian[outerLupe][innerLupe] = stiff[outerLupe][innerLupe];
        }
    for(outerLupe=0;outerLupe<=N;++outerLupe)
    {
        jacobian[outerLupe][outerLupe] -= gaussWeights[outerLupe];
        baseFunc[outerLupe] = 1.0*gaussWeights[outerLupe];
    }

}
