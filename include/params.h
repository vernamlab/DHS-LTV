#ifndef PARAMS_H
#define PARAMS_H

#include "general.h"

class Params
{
    public:
        Params(){delta_ = DELTA; n_ = N; logq_ = SIGMA; hermite_ = 100;};
        virtual ~Params();

        double  ComputeHermite(){hermite_ = pow(2, (logq_-4)/(4*n_)); return hermite_;};
        bool    CheckSecurity(){return (hermite_ < delta_);};
        bool    CheckNoiseGrowth();
    protected:
    private:
            double delta_;
            double hermite_;
            int n_;
            int logq_;
};

#endif // PARAMS_H
