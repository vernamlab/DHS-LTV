#include "fft_mult.h"

long GetPrimeBound(int deg_a, int deg_b, int bitsize_a, int bitsize_b, int num_of_additions){
	long res;
	int log_add;

	log_add = NumBits(num_of_additions)+1;
	res 	= 2 + NumBits(min(deg_a, deg_b)+1) + bitsize_a + bitsize_b + log_add;
	return res;
}

long GetPrimeNumber(long bound, ZZ &prod){

	long nprimes;
	zz_pBak bak;
	bak.save();
	for (nprimes = 0; NumBits(prod) <= bound; nprimes++) {
		UseFFTPrime(nprimes);
		      mul(prod, prod, GetFFTPrime(nprimes));
   }
	bak.restore();
	return nprimes;
}

void CalculateFFTValues(fftRep *R1, ZZX &a, int &prime_num, int deg_b){

	zz_pBak bak;
	bak.save();

	zz_pX A;
	long k, d;

	d = 2*deg_b;
	k = NextPowerOfTwo(d+1);
	for(int i=0; i<prime_num; i++){
		zz_p::FFTInit(i);
		conv(A, a);
	    TofftRep(R1[i], A, k);
	}

	bak.restore();
}

void fftRepClear(fftRep *R, int &priNum){
	long* p;
	for(int i=0; i<priNum; i++){
		p = &(R[i].tbl[0][0]);
		for(int j=0; j<(1L<<R[i].k); j++){
			p[j] = 0;
		}
	}
}

void fftRepSetSize(fftRep *R, int &priNum, long &fftSize){
	for(int i=0; i<priNum; i++){
		zz_p::FFTInit(i);
		R[i].SetSize(fftSize);
	}
}

void mymult(){
	ZZX mya, myb, c0, c1, x;
	ZZ q;
	int k = to_long(euler_toient(to_ZZ(Modulus_M)));

	GenPrime(q, Max_Prime);
	RandomPolyGen(mya, k, 1, q);
	RandomPolyGen(myb, k, Max_Prime, q);

	long da = deg(mya);
	long db = deg(myb);
	long bound = 2 + NumBits(min(da, db)+1) + MaxBits(mya) + MaxBits(myb);


	ZZ prod;
	set(prod);
	int prime_num = GetPrimeNumber(bound, prod);
	cout << prime_num << endl;


	long mk = NextPowerOfTwo(2*da+1);

	zz_p::FFTInit(0);
	long p = zz_p::modulus();

	fftRep R1[prime_num];
	fftRep R2[prime_num];
	fftRep R3[prime_num];
	fftRep R4[prime_num];

	int size = 256;
	fftRep Rm[prime_num][size];

	for(int i=0; i<prime_num; i++)
		for(int j=0; j<size; j++)
			Rm[i][j].SetSize(mk);

	for(int i=0; i<prime_num; i++){
		zz_p::FFTInit(i);
		R1[i].SetSize(mk);
		R2[i].SetSize(mk);
		R3[i].SetSize(mk);
		R4[i].SetSize(mk);
	}

	myTimer tm;
	tm.Start();
	CalculateFFTValues(R1, mya, prime_num, db);
	tm.Stop();
	tm.ShowTime("My FFT:\t");

	CalculateFFTValues(R2, myb, prime_num, db);

	tm.Start();


	for(int i=0; i<prime_num; i++)
		for(int j=0; j<size; j++)
			Rm[i][j] = R2[i];

	for(int j=0; j<size; j++){
		CalculateFFTValues(R1, mya, prime_num, db);
		for(int i=0; i<prime_num; i++){
			zz_p::FFTInit(i);
			mul(R3[i], R1[i], Rm[i][j]);
			add(R4[i], R4[i], R3[i]);
		}
	}
	CalculateFFTValues(R4, myb, prime_num, db);

	tm.Stop();
	tm.ShowTime("My FFT:\t");
}


