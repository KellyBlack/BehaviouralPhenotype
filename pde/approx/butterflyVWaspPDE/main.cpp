#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <thread>

#include "numericaltrials.h"

//#define OUTPUTFILE "approximation.csv"
#define BINARYOUTPUTFILE "approximation"
#define SKIP_PRINT_UPDATE 10000
#define SKIP_FILE_SAVE 150
#define NUMBER_TIME_LOOP 3000000
#define MAX_NEWTON_STEPS 50
#define LEGENDRE_POLY_DEGREE 30
#define MAX_DELTA_NORM 0.0001

int saveSeparateRunsByM(double mu,double c,double g,double d,
                        double lowM,double highM,double stepM,double dt);

int main()
{
    // Set up the temporal variables.
    double dt = 0.0001;

    //NumericalTrials *trials = new NumericalTrials();

    // Define the default values of the parameters.
    double mu = 0.01;
    double c  = 0.7;
    double g  = 2.0;
    double d  = 1.0;
    //double m  = 0.2;

    std::cout << "Starting" << std::endl;
    saveSeparateRunsByM(mu,c,g,d,0.2,1.3,0.2,dt);
    //saveSeparateRunsByM(mu,c,g,d,0.4,0.85,0.2,dt);
    std::cout << "Done" << std::endl;

    //delete trials;
    return(0);
}


#define NUMBER_THREADS 10
int saveSeparateRunsByM(double mu, double c, double g, double d,
                        double lowM, double highM, double stepM, double dt)
{

    int num = 0;
    double m = lowM;
    std::vector<std::thread> processes;
    std::vector<NumericalTrials*> trials;
    while(m<=highM)
    {
        processes.clear();
        trials.clear();
        num = 0;
        while((m<=highM)&&(num<NUMBER_THREADS))
        {
            std::ostringstream filename("");
            filename << BINARYOUTPUTFILE
                     << "-m-" <<  std::setw(6) << std::fixed << std::setprecision(4) << std::setfill('0') << m
                     << "-mu-" << mu
                     << ".bin";
            std::cout << "Writing to " << filename.str() << std::endl;

            NumericalTrials *trial = new NumericalTrials();
            trials.push_back(trial);
            processes.push_back(
                        std::thread(&NumericalTrials::approximateSystem,trial,
                                    mu,c,g,d,m,
                                    dt,NUMBER_TIME_LOOP,
                                    LEGENDRE_POLY_DEGREE,
                                    MAX_DELTA_NORM,MAX_NEWTON_STEPS,
                                    filename.str(),
                                    SKIP_PRINT_UPDATE,
                                    SKIP_FILE_SAVE)
                        );

            std::cout << "Moving along " << m << std::endl;
            m += stepM;
            num += 1;
        }

        std::vector<std::thread>::iterator eachProcess;
        for(eachProcess=processes.begin();eachProcess!=processes.end();++eachProcess)
        {
            eachProcess->join();
        }

        std::vector<NumericalTrials*>::iterator eachTrial;
        for(eachTrial=trials.begin();eachTrial!=trials.end();++eachTrial)
        {
            delete *eachTrial;
        }
        //th.join();
    }

    return(1);
}
