#ifndef LTV_H_
#define LTV_H_

#include <fstream>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pE.h>

#include "fft_mult.h"

using namespace std;
NTL_CLIENT

struct eval_key{
	ZZX *key;
};

class ltv{
public:

///////////////////////////////////////////////////////////////////////////////
	ltv(int dm, int num, int max_bit, int tableS, int wordSize);
	ltv(ZZ my_p, ZZ my_B, int my_N, int my_dm, int num, int max_bit, int tableS, int wordSize);

////// Primitive Functins
		ZZX 	Prim_Encrypt(ZZX m, int index);
		ZZX		Prim_Decrypt(ZZX c, int index);
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
		void 	Arith_PolyCoeffReduce(ZZX &out, ZZX &in, int index);
		ZZX		Arith_MulModZZX(ZZX &a, ZZX &b, int index);
		ZZX		Arith_AddModZZX(ZZX &a, ZZX &b, int index);
///////	Pointers
		ZZ		*Pointer_Q(int i);
		ZZX		*Pointer_PolyModulus();
		ZZX		*Pointer_PublicKey();
		ZZX		*Pointer_SecretKey(int i);
		ZZX		*Pointer_EvalKey(int i);
		GlobalParam *Pointer_GP();

private:
	int word;
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
	int rand_seed, message_rand, max_bitsize, word_blocksize; //

/////// Private Functions
	void	GetPolyBitsIndex(ZZX &bits, ZZX &cip, int &bitIndex);								// Get the selected index bits of the coefficients of the polynomials
	void	GetPolyWordsIndex(ZZX &bits, ZZX &cip, int &bitIndex);								// Get the selected index words of the coefficients of the polynomials
	void 	FFTmult(fftRep *res, fftRep *R1, fftRep *R2, int &priNum);							// FFT multiplication
	void 	FFTadd(fftRep *res, fftRep *R1, fftRep *R2, int &priNum);							// FFT addition
	void 	FromCrt2Poly(ZZX &res, zz_pX *polys, int &priNum);									// Convert from CRT rep. to Polynomial
	void	ConvertKeystoFFT(eval_key &ek2, FFTPolyList &fftkey, int &bitsize, int &num_add);	// Convert evaluation keys to FFT representation
	void	compute_eval_onerelin(eval_key &ek2, ZZX tppk, ZZX tpek, int index);				// Compute Relin
	ZZX 	ComputeFastCycModulus(int n);														// Compute the cyclotomic polynomial for given M
	void 	SetFastPolyModulus();																// Set polynomials for fast polynomial modular reduction
	int 	NumofWords(int bs);																	// Compute bitsize/wordsize + 1
};
///////////////////////////////////////////////////////////////////////////////
////// Functions used by ltv member functions
void 	clear(ZZX &x, int size);
void 	find_inverse(ZZX &f_inv, ZZX &f, ZZ &q, int &degree, bool &isfound, ZZX modulus);
void 	coeff_reduction(ZZX &out, ZZX &in, ZZ &q, int &degree);
void 	getClosestMod2(ZZ& out, ZZ& in, ZZ& p1, ZZ& p2);
int 	MobuisFunction(int n);
///////////////////////////////////////////////////////////////////////////////

#endif /* LTV_H_ */

