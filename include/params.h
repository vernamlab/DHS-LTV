#ifndef PARAMS_H
#define PARAMS_H

#include "general.h"

class Params
{
    public:
        Params(double delta, int degree, int qsize);
        virtual ~Params(){};

        double  ComputeHermite();
        bool    CheckSecurity();
        bool    CheckNoiseGrowth();
    protected:
    private:
            double delta_;
            double hermite_;
            int n_;
            int logq_;
};

#endif // PARAMS_H
