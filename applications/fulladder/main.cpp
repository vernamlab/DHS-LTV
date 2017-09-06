#include "ltv.h"

using namespace std;

void SetMessage(ZZX *bits, uint8_t &m, int size){
	for (int i = 0; i < size; ++i){
		bits[i] = (m >> i) % 2;
	}
}

void EncMessage(ZZX *bits, ltv &n, int size){
	for (int i = 0; i < 4; ++i){
		bits[i] = n.Prim_Encrypt(bits[i], 0);
	}
}

void FullAdder(ZZX *Sum, ZZX *CarryOut, ZZX *Cin, ZZX *Bits_a, ZZX *Bits_b, ltv &n, int size){
	ZZX Temp;
	Sum[0] = Bits_a[0] + Bits_b[0];
	Temp = Bits_a[0] * Bits_b[0];
	Temp = n.Arith_PolyCoeffReduce(Temp,0);
	Temp = n.Prim_RelinRingFFT(Temp, 0);
	Temp = n.Arith_PolyCoeffReduce(Temp,0);
	Temp = n.Prim_ModSwitch(Temp, 0);
	CarryOut[0] = Sum[0] * Cin[0];
	CarryOut[0] = n.Arith_PolyCoeffReduce(CarryOut[0],0);
	CarryOut[0] = n.Prim_RelinRingFFT(CarryOut[0], 0);
	CarryOut[0] = n.Arith_PolyCoeffReduce(CarryOut[0],0);
	CarryOut[0] = n.Prim_ModSwitch(CarryOut[0], 0);
	CarryOut[0] += Temp;
	Sum[0] += Cin[0];
}

int main(){

	int M 			= Modulus_M;			// Mth Cyclotomic Polynomial, n=euler_phi(M)
	ZZ p 			= to_ZZ("2");	// Message modulus
	int wordsize	= 16;			// radixsize. when computing relinearization we attach wordsize bits. larger makes it faster, smaller evaluaation keysize. however cause larger noise. usually 16 is limit, after that need to cut more bits.
	int numLevel	= Num_Primes-1;			// Number of levels
	int numPrimes	= numLevel+1;	// Num of Primes is equal to Number of levels + 1
	int logq		= Dif_Prime;			// Prime size which is also noise cutting size.
	int Maxlogq		= numPrimes*logq; //Maximum bitsize of the modulus q.
	int tableNum 	= 1;

	ltv n(M, numPrimes, Maxlogq, tableNum, p, wordsize);
	n.Func_SetModulus();
	n.Func_ModulusFindRing(numPrimes, Maxlogq, logq, *n.Pointer_PolyModulus());
	n.Func_ComputeKeysRingRelinFFT(numPrimes, tableNum);

	uint8_t CarryIn = 0x00;
	uint8_t a = 0x03, b = 0x06;
	ZZX Bits_a[4], Bits_b[4], Sum[4], CarryOut[4],  Cin[1];
	SetMessage(Bits_a, a, 4);
	SetMessage(Bits_b, b, 4);	
	SetMessage(Cin, CarryIn, 1);

	EncMessage(Bits_a, n, 4);
	EncMessage(Bits_b, n, 4);
	EncMessage(Cin, n, 1);
	
	FullAdder(Sum, CarryOut, Cin, Bits_a, Bits_b, n, 4);

	CarryOut[0] = n.Arith_PolyCoeffReduce(CarryOut[0], 1);
	CarryOut[0] = n.Prim_Decrypt(CarryOut[0], 1);
	Sum[0] = n.Prim_Decrypt(Sum[0], 0);

	cout << CarryOut[0] << Sum[0] << endl;


	return 0;
}



