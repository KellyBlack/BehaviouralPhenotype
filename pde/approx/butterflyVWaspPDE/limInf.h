#ifndef LIMINF_H
#define LIMINF_H

#include <fstream>
#include <sstream>
#include <iostream>

template <class number>
class LimInf
{

public:
    LimInf(number initial,bool maximum=true);
    number extreme();
    void setExtreme(number val);
    void operator=(number val);
private:

    bool testDirection(number currentVal,number prevVal);

    bool searchMax;
    bool trending;
    long iterations;
    number extremeValue;
    number prevExtremeValue;
    number prevValue;
};

template <class number>
LimInf<number>::LimInf(number initial,bool maximum)
{
    searchMax = maximum;
    trending = false;
    iterations = 0;
    extremeValue = initial;
    prevValue = initial;
    prevExtremeValue = initial;
}

template  <class number>
void LimInf<number>::setExtreme(number val)
{
    extremeValue = val;
    prevExtremeValue = val;
}

template <class number>
void LimInf<number>::operator=(number val)
{
    //std::cout << "Checking " << val << ", " << prevExtremeValue << ": ";
    //if(searchMax)
    {
        // We are trying to keep track of the largest value in the cycle.
        if(trending)
        {
            // The value of the function is assumed to be increasing.
            if(testDirection(val,prevValue))
            {
                // The value of the function is getting bigger.
                extremeValue = val;
                iterations = 0;             // This the last time the function was increasing.
            }

            // The function is decreasing. Is it a minor blip or part of a
            // long term trend?
            else if(++iterations > 100)
            {
                // It has been decreasing for a while. Assume it is part of a long term trend,
                // and the function is decreasing.
                prevExtremeValue = extremeValue;
                trending = false;
                iterations = 0;
            }
            //std::cout << "Increasing  " << extremeValue << std::endl;

        }

        else
        {
            // The value of the function is assumed to be decreasing in the big picture.
            //std::cout << "Decreasing  " << extremeValue << std::endl;
            if(!testDirection(val,prevValue))
            {
                // The function is decreasing.
                iterations = 0;
            }
            else if(++iterations > 100)
            {
                // It has been increasing for a while. Assume it is part of a long term trend,
                // and the function is now increasing
                trending = true;
                iterations = 0;
                extremeValue = val;
            }

        }

        prevValue = val; // keep track of the previous value.

    }

}


template <class number>
bool LimInf<number>::testDirection(number currentVal,number prevVal)
{
    if(searchMax)
    {
        return(currentVal>prevVal);
    }
    return(currentVal<prevVal);
}

template <class number>
number LimInf<number>::extreme()
{
    return(prevExtremeValue);
}


#endif // LIMINF_H
