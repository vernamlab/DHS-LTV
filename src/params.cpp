#include "params.h"

Params::Params(double delta, int degree, int qsize)
{
    delta_ = delta;
    n_ = degree;
    logq_ = qsize;
    ComputeHermite();
}

double Params::ComputeHermite()
{
    hermite_ = pow(2, (logq_-4)/(4*n_));
    return hermite_;
}

bool Params::CheckSecurity()
{
    return (hermite_ < delta_);
}

bool Params::CheckNoiseGrowth()
{

}
