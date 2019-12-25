#include "approximationbase.h"

ApproximationBase::ApproximationBase()
{

}

ApproximationBase::MessageInformation* ApproximationBase::findReturnedProcessParameters(
        ApproximationBase::MaxMinBuffer msgValue,
        std::vector<ApproximationBase::MessageInformation*> processes)
{

    std::vector<ApproximationBase::MessageInformation*>::iterator eachProcess;

    // need to find this thread and join it.
    for(eachProcess=processes.begin();(eachProcess!=processes.end());++eachProcess)
    {
        if((*eachProcess)->which == msgValue.which)
        {
            // This is the thread from which this process sprang forth.
            // record the values that were passed and clean up the mess.
            (*eachProcess)->mu    = msgValue.mu;
            (*eachProcess)->theta = msgValue.theta;
            (*eachProcess)->c     = msgValue.c;
            (*eachProcess)->g     = msgValue.g;
            (*eachProcess)->d     = msgValue.d;
            (*eachProcess)->m     = msgValue.m;
            (*eachProcess)->time  = msgValue.endTime;
            (*eachProcess)->maxButterfly = msgValue.maxButterfly;
            (*eachProcess)->minButterfly = msgValue.minButterfly;
            (*eachProcess)->maxWasp      = msgValue.maxWasp;
            (*eachProcess)->minWasp      = msgValue.minWasp;
            return(*eachProcess);
        }
    }

    return(nullptr);

}

