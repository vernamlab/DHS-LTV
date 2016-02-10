#ifndef GENERAL_H_
#define GENERAL_H_

#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZXFactoring.h>
#include <time.h>
#include "def.h"

using namespace std;
NTL_CLIENT

//////////////////////////////////////////////////////////////////////////////////
typedef enum SetupType
{
    automatic,
    manual
}SetupType;

typedef enum RingType
{
    cyclotomic,
    xn_1
}RingType;

typedef enum ReductionType
{
    fast,
    barrett,
    ntl
}ReductionType;

typedef enum HomType
{
    fhe,
    lhe,
    she
}HomType;

typedef enum Flag
{
    on,
    off
}Flag;

typedef struct Ring
{
    vec_ZZ q_;
    ZZX poly_;
}Ring;


typedef struct FFTReps
{
    int block_cnt_ = 0;             // For each ek (block_cnt = #relin blocks or #ek's)
    int rep_cnt_ = 0;               // For each fft_prime (rep_cnt = #fft primes)
    fftRep **rep_list_ = NULL;      // [block_cnt * rep_cnt] matrix
}FFTReps;


typedef struct FFTEvalKey
{
    int length_ = 0;            // For each level
    FFTReps *level_ = NULL;
}FFTEvalKey;


///////////////////////////////////////////
// FREE FUNCTIONS ~ Polynomial Operations
void 	Clear(ZZX &x);
void    RandPolyBalanced(ZZX &out, int n, int q);
void    RandPolyBalanced(ZZX &out, int n, const ZZ &q);
void 	RandPolyBounded(ZZX &x, int n, int q);
void 	RandPolyBounded(ZZX &x, int n, const ZZ &q);
void    CyclotomicPoly(ZZX &out, int m);
void    XN_1(ZZX &out, int n);
void 	FindInverse(ZZX &out, const ZZX &in, const ZZX& poly_mod, const ZZ &coeff_mod, bool &found);
void    BarrettReduction(ZZX &out, const ZZX &in, const ZZX &mod, const ZZX &u);
void	BarrettReduction(ZZX &out, const ZZX &in, const ZZX &mod);
void    Div_high_low(ZZX &high, ZZX &low, ZZX &t, int low_deg, int high_deg);
void    fastReduction(ZZX &out, const ZZX &in, const ZZX &mod, const vec_ZZX &helpers);
void    fastReduction(ZZX &out, const ZZX &in, const ZZX &mod);
void	PolyReduce(ZZX &out, const ZZX &in, RingType type, const ZZX &mod);
void	CoeffReduce(ZZX &out, const ZZX &in, const ZZ &mod);
void	CoeffReduce(ZZX &out, const ZZX &in, int q);
void    CoeffBalance(ZZX &out, const ZZX &in, const ZZ &q);
void    CoeffDivide(ZZX &out, const ZZX &in, const ZZ &q);
void    CoeffParityMatch(ZZX &out, const ZZX &in, int q);
void 	PolyCoeffReduce(ZZX &out, const ZZX &in, RingType type, const ZZX &poly_mod, const ZZ &coeff_mod);

bool 	FindInverse(ZZX &out, const ZZX &in, const ZZX& poly_mod, const ZZ &coeff_mod);
void    ComputeFFT(fftRep *rep_list, const ZZX &poly, int prime_cnt, int deg);
void    FFTRepClear(fftRep *rep_list, int length);
void    FFTRepSetK(fftRep *rep_list, int length, int k);
void    SeparateBlocks(vec_ZZX &out, const ZZX &in, long blocksize);
void    FFTtoPoly(ZZX &out, vec_zz_pX &in, int rep_length, int n);

// Big Number Operations
long    GetPrimeNumber(long bound, ZZ &prod);
void    FindPrimeCongOne(ZZ &out, const ZZ &mod);
void    FactorizeNumber(vec_ZZ &factors, vec_ZZ &powers, const ZZ &in);
void 	FindGenerator(ZZ &out, const vec_ZZ &factors, const ZZ &p);
void 	FindNthRootOfUnity(ZZ &out, const ZZ &g, const ZZ &p, const ZZ &n);
void 	getClosestMod2(ZZ &out, const ZZ &in, const ZZ &p1, const ZZ &p2);
void 	getClosestModP(ZZ &out, const ZZ &in, const ZZ &p1, const ZZ &p2, const ZZ &modp);
void    EulerToient(ZZ &out, const ZZ &in);

int 	MobuisFunction(int n);
int     FindPrimeCongOne(int mod);
int     EulerToient(int in);
int		ComputeFactorDegree(int m);



//////////////////////////////////////////////////////////////////////////////////

class myTimer{
public:
	void Start();
	void Stop();
	double GetTime();
	void ShowTime(string s);
private:
	clock_t startTime;
	clock_t stopTime;
};


#endif

















