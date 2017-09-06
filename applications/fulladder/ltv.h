#ifndef ltv_H_
#define ltv_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZXFactoring.h>

#include "general.h"
#include "fft_mult.h"
#include <stdlib.h>
#include <ctime>

using namespace std;
NTL_CLIENT

struct eval_key{
	ZZX *key;
};

class ltv{
public:

///////////////////////////////////////////////////////////////////////////////
	ltv(int dm, int num, int max_bit, int tableS, ZZ qq,  int wordSize);
	ltv(ZZ my_p, ZZ my_B, int my_N, int my_dm, int num, int max_bit, int tableS,  int wordSize);

////// Primitive Functins
		ZZX 	Prim_Encrypt(ZZX m, int index);
		ZZX		Prim_Decrypt(ZZX c, int index);
		ZZX		Prim_Decrypt2(ZZX c, int index);
		ZZX 	Prim_ModSwitch(ZZX &x, int index);
		ZZX 	Prim_ModSwitch_LongPoly(ZZX &x, int size, int index);
		ZZX 	Prim_RelinRingFFT(ZZX &c, int index); 					// Relinzation
		ZZX 	Prim_RelinRingFFT(ZZX &c, int index, int tableIndx);	// Relinzation for the set index
		void 	Prim_ReduceKeysLevel(int level);						// Update Table to next level
		void 	Prim_ReduceKeysLevel(int level, int tblIndex);			// Update Table to next level with index
////// Some Useful Functions
		ZZX 	Func_Sample();													// Create Sample poly with B-bound
		ZZX 	Func_CreateMessage(int size);									// Create Random Message
		void 	Func_SetModulus();												// Set Modulus
		void   	Func_ModulusFindRing(int num, int max_bit, int diff, ZZX modu);	// Compute modulus
		void	Func_ComputeKeysRingRelinFFT(int num, int tblsize);				// Compute Eval keys in FFT form
//////	Arithmetic Operations
		void	Arith_PolyReduce(ZZX &out, ZZX &in);
		void	Arith_CoeffReduce(ZZX &out, ZZX &in, int index);
		ZZX 	Arith_PolyCoeffReduce(ZZX &in, int index);
		ZZX		Arith_MulModZZX(ZZX &a, ZZX &b, int index);
		ZZX		Arith_AddModZZX(ZZX &a, ZZX &b, int index);
///////	Pointers
		ZZ		*Pointer_Q(int i);
		ZZX		*Pointer_PolyModulus();
		ZZX		*Pointer_PublicKey();
		ZZX		*Pointer_SecretKey(int i);
		ZZX		*Pointer_EvalKey(int i);
		GlobalParam *Pointer_GP();
//////// IO Handle
		void 	IO_ReadModulus(fstream &fs);
		void 	IO_ReadQs(fstream &fs, int s);
		void	IO_ReadKeys(fstream &fs, int s);
		ZZX key_f, key_fi;
private:
	int word, logq;
	int N, degree_m;  		// poly degree, euler phi_m
	ZZ 	p, q, B;			// message modulus, encryption modulus, noise bound
	ZZX modulus;			// polynomial modulus


	ZZ  	*q_list;		// list of q modulus for each level
	ZZX 	*pk, *sk;		// list of public and secret keys for each level

	struct eval_key *ek2;	// list of evaluation keys for some table size
	FFTPolyList *fftKeys;	// list of eval keys in fft form
	int tableSize;			// table size

	GlobalParam gp;			// Global param struct to store some parameters
	myReduction myr;		// class to complete fast polynomial reduction
	bool reducsetflag;		//
	int rand_seed, message_rand, max_bitsize, init_bitsize, word_blocksize; //

/////// Private Functions
	void	GetPolyBitsIndex(ZZX &bits, ZZX &cip, int &bitIndex);
	void	GetPolyWordsIndex(ZZX &bits, ZZX &cip, int &bitIndex);
	void 	FFTmult(fftRep *res, fftRep *R1, fftRep *R2, int &priNum);
	void 	FFTadd(fftRep *res, fftRep *R1, fftRep *R2, int &priNum);
	void 	FromCrt2Poly(ZZX &res, zz_pX *polys, int &priNum);
	void	ConvertKeystoFFT(eval_key &ek2, FFTPolyList &fftkey, int &bitsize, int &num_add);
	void	compute_eval_onerelin(eval_key &ek2, ZZX tppk, ZZX tpek, int index);
	ZZX 	ComputeFastCycModulus(int n);
	void 	SetFastPolyModulus();
	int 	NumofWords(int bs);
};
///////////////////////////////////////////////////////////////////////////////
////// Functions used by ltv member functions
void 	clear(ZZX &x, int size);
void 	find_inverse(ZZX &f_inv, ZZX &f, ZZ &q, int &degree, bool &isfound, ZZX modulus);
void 	coeff_reduction(ZZX &out, ZZX &in, ZZ &q, int &degree);
void 	coeff_reduction_q_2(ZZX &out, ZZX &in, ZZ &q, int N);
void 	getClosestMod2(ZZ& out, ZZ& in, ZZ& p1, ZZ& p2);
void 	RoundModP(ZZ& out, ZZ& in, ZZ& p1, ZZ& p2, ZZ& modp);
ZZ 		ComputePR(ZZ q, ZZ m);
vec_ZZ 	findFactors(ZZ f);
int 	MobuisFunction(int n);
///////////////////////////////////////////////////////////////////////////////

#endif /* ltv_H_ */





