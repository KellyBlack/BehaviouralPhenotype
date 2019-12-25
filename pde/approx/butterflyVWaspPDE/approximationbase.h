#ifndef APPROXIMATIONBASE_H
#define APPROXIMATIONBASE_H

#include <thread>
#include <vector>


class ApproximationBase
{
public:
    ApproximationBase();

protected:

    // Set up the data structure that will be used to keep track of the processes
    // and returned information associated with each process.
    struct MessageInformation
    {
        long which;
        double theta;
        double mu;
        double c;
        double g;
        double d;
        double m;
        double time;
        double maxButterfly;
        double minButterfly;
        double maxWasp;
        double minWasp;
        std::thread *process;
        ApproximationBase *trial;
    };

    struct MaxMinBuffer
    {
      long mtype;
      long which;
      double theta;
      double c;
      double g;
      double d;
      double m;
      double mu;
      double endTime;
      double maxButterfly;
      double minButterfly;
      double maxWasp;
      double minWasp;
    };

    MessageInformation* findReturnedProcessParameters(
            MaxMinBuffer msgValue,
            std::vector<MessageInformation*> processes);

};

#endif // APPROXIMATIONBASE_H
