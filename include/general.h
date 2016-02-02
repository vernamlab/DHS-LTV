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
typedef enum RingType
{
    cyclotomic  = 1,
    xn_1        = 2
}RingType;

typedef enum ReductionType
{
    barrett     = 1,
    yarkin      = 2,
    ntl         = 3
}ReductionType;

typedef enum HomType
{
    fhe     = 1,
    she     = 2
}HomType;

typedef enum BatchFlag
{
    on      = 1,
    off     = 2
}BatchFlag;

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



/*
typedef struct CipherText
{
    ZZX ct_;
    bool fresh_;
    int level_;
    int multcnt_;
}CipherText;
*/


class CipherText
{
public:
            CipherText(){set_d(to_ZZ(1)); set_ct(0); set_fresh(true); set_level(0); set_multcnt(0);}
            // Operator Overloading
            CipherText& operator=(const CipherText &rhs);
            CipherText& operator+=(const CipherText &rhs);
            const CipherText operator+(const CipherText& rhs) const;
            CipherText& operator-=(const CipherText &rhs);
            const CipherText operator-(const CipherText& rhs) const;

            CipherText& operator+=(const ZZ &rhs);
            const CipherText operator+(const ZZ& rhs) const;
            CipherText& operator*=(const ZZ &rhs);
            const CipherText operator*(const ZZ& rhs) const;

            friend std::ostream& operator<< (std::ostream&, const CipherText& ct);

            const ZZX&  ct() const{ return ct_; };
            const ZZ&   d() const{ return d_; };
            bool        fresh()const{return fresh_;};
            int         level()const{return level_;};
            int         multcnt()const{return multcnt_;};

            void    set_ct(int i){ct_ = i;};
            void    set_ct(const ZZX &ct){ct_ = ct;};
            void    set_fresh(bool f){fresh_ = f;};
            void    set_level(int l){level_ = l;};
            void    set_multcnt(int cnt){multcnt_ = cnt;};
            void    set_d(const ZZ &d){d_ = d;};

            void    coeff_reduce(const ZZ &mod);
            void    poly_reduce(RingType type, const ZZX &mod);

private:
            ZZX ct_;
            bool fresh_;
            int level_;
            int multcnt_;
            ZZ  d_;
};


class GlobalParams
{
public:
        GlobalParams();
        GlobalParams(RingType type, BatchFlag flag, int p, int degree);
        //~GlobalParams();

        //void        Init();
        void            Init(RingType type, BatchFlag flag, int p, int degree);


        HomType         hom_type()const{return hom_type_;};
        RingType        ring_type()const{return ring_type_;};
        BatchFlag       flag()const{return batch_flag_;};
        ReductionType   reduc_type()const{return reduc_type_;};
        int             p()const{return p_;};
        int             n()const{return n_;};
        int             m()const{return m_;};
        int             k()const{return k_;};
        int             i()const{return i_;};
        int             b()const{return b_;};
        int             l()const{return l_;};
        int             d()const{return d_;};
        int             r()const{return r_;};
        int             w()const{return w_;};



        void        set_hom_type(HomType type){hom_type_ = type;};
        void        set_ring_type(RingType type){ring_type_ = type;};
        void        set_batch_flag(BatchFlag flag){batch_flag_ = flag;};
        void        set_reduc_type(ReductionType type){reduc_type_ = type;};
        void        set_p(int p){p_ = p;};
        void        set_n(int n){n_ = n;};
        void        set_m(int m){m_ = m;};
        void        set_k(int k){k_ = k;};
        void        set_i(int i){i_ = i;};
        void        set_l(int l){l_ = l;};
        void        set_r(int r){r_ = r;};
        void        set_b(int b){b_ = b;};
        void        set_d(int d){d_ = d;};
        void        set_w(int w){w_ = w;};


private:
        HomType         hom_type_;
        RingType        ring_type_;
        BatchFlag       batch_flag_;
        ReductionType   reduc_type_;
        int             p_, n_, m_, k_, i_, b_, l_, d_, r_, w_;


};

/*
typedef struct{
    // Application Params
    //int L;      // user input bit cnt
    //int R;      // output cnt

    RingType t;

    // Security Params
    int k;      // cut bit-size
    int l;      // how many times we will cut?
    int i;      // the last q size
    //int maxsize;
	//ZZ* q_list; // q[i] = q[l] + k^i
	int p;      // plaintext space, 2 or 257?
	int B;       // noise bound, usually 1

	// Cyclotomic Poly Params
	int m;  // cyc m
	int n;  // cyc degree n = Phi(m)
	int d;  // fact degree
	int r;  // fact cnt n = d.r, slot count
}GlobalParam;

void 	GlobalParamInit(GlobalParam *gp);
//void    GlobalParamKill(GlobalParam &gp);
*/

// FHE Related
//void    FindQList(vec_ZZ &q_list, const GlobalParam *gp);

///////////////////////////////////////////
// Polynomial Operations
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
void    YarkinReduction(ZZX &out, const ZZX &in, const ZZX &mod, const vec_ZZX &helpers);
void    YarkinReduction(ZZX &out, const ZZX &in, const ZZX &mod);
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

















