#ifndef BUTTERFLIES_H
#define BUTTERFLIES_H

#include "pdeSolver.h"

class Butterflies : public PDESolver
{
public:
    Butterflies(int number=0);

    virtual void buildJacobian();

protected:
    ~Butterflies();

private:

};

#endif // BUTTERFLIES_H
