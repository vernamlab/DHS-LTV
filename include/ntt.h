#ifndef NTT_H_
#define NTT_H_

#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include "general.h"

using namespace std;
using namespace NTL;

class BatcherNTT{

public:
            //BatcherNTT(GlobalParam *gp);
            BatcherNTT(int p, int degree);

            void    Init(int p, int degree);
            ZZ      p(){return p_;};
            void    FindP(int p);


            void    CooleyTukey(ZZX &out, const ZZX &in, const ZZ &r, int l);
            void  	Encode(ZZX &out, const ZZX &in);
            void  	Decode(ZZX &out, const ZZX &in);
            void 	fft_core(ZZX &out, const ZZX &in, const ZZ &r, int l);


            void    PrintFactors();
            void    PrintP(){cout << "p : " << p_ << endl;};
            void    PrintDegree(){cout << "n : " << n_ << endl;};
            void    PrintGenerator(){cout << "g : " << g_ << endl;};
            void    PrintNthRootOfUnity(){cout << "alpha : " << alpha_ << endl;};
            void    PrintInverseNthRootOfUnity(){cout << "alpha^{-1} : " << inv_alpha_ << endl;};


private:
            int         n_;
            ZZ          inv_n_, p_;
            ZZ          g_, alpha_, inv_alpha_;     // g_:generator mod p,
                                                    // alpha_: nth primitive root of unity mod p
            vec_ZZ      factors_, powers_;
};
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

class BatcherCRT{
public:

    //BatcherCRT(const GlobalParam *gp);
    BatcherCRT(int cyc_m, int deg, int prime_p, int fac_deg, int fac_cnt);
    BatcherCRT(const ZZ &cyc_m, const ZZ &deg, const ZZ &prime_p, int fac_deg, int fac_cnt);

	void	Encode(ZZX &out, const ZZX &in);
	void    Decode(ZZX &out, const ZZX &in);

private:
	string     mod_file;
	string     factor_file;
	string     product_file;
	int         d;
	int 		r;
	int         m;
	int         n;
	int         p;
};

#endif /* NTT_H_ */
