#ifndef FFT_MULT_H_
#define FFT_MULT_H_

#include "def.h"
#include "general.h"
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/RR.h>
#include <NTL/lzz_pX.h>
#include <NTL/FFT.h>

NTL_CLIENT

///////////////////////////////////////////////////////////////////////////////////////////////
typedef struct{
	int 	rep_num;
	fftRep *poly_rep;
}PolyFFTRep;

typedef struct{
	int			 polynum;
	PolyFFTRep 	*polyList;
}FFTPolyList;
///////////////////////////////////////////////////////////////////////////////////////////////

long GetPrimeBound(int deg_a, int deg_b, int bitsize_a, int bitsize_b, int num_of_additions);
long GetPrimeNumber(long bound, ZZ &prod);
void CalculateFFTValues(fftRep *R1, ZZX &a, int &prime_num, int deg_b);
void fftRepClear(fftRep *R, int &priNum);
void fftRepSetSize(fftRep *R, int &priNum, long &fftSize);
////
void mymult();



#endif
