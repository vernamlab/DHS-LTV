#ifndef CRT_H_
#define CRT_H_

#include <fstream>
#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pXFactoring.h>

using namespace std;
NTL_CLIENT

class myCRT{
public:
	void SetModulus(ZZX &m);
	void ComputeFactors(int f_degree, int f_size);
	void CalculateMs();
	void CalculateNs();
	void CalculateMxNs();
	bool CRTTestNMs();
	ZZX	 EncodeMessageMxN(ZZX &mess);
	ZZX	 DecodeMessage(ZZX &mess);

	ZZX 		num2ZZX(int num);
	vec_ZZ_pX 	ReturnM();
	vec_ZZ_pX 	ReturnN();
	vec_ZZ_pX 	ReturnMxN();
	vec_ZZ_pX 	Returnfactors();
	ZZ_pX	  	Returnmodulus();
	int 	  	Returnsize();

	void 		IOReadAll(fstream &fs);
	void 		IOReadModulus(fstream &fs);
	void 		IOReadSize(fstream &fs);
	void 		IOReadFactors(fstream &fs);
	void 		IOReadMs(fstream &fs);
	void 		IOReadNs(fstream &fs);
	void 		IOReadMxNs(fstream &fs);

	void 		IOReadAll();
	void 		IOReadModulus();
	void 		IOReadSize();
	void 		IOReadFactors();
	void 		IOReadMs();
	void 		IOReadNs();
	void 		IOReadMxNs();

private:
	vec_ZZ_pX 	M;
	vec_ZZ_pX 	N;
	vec_ZZ_pX	MxN;
	vec_ZZ_pX 	factors;
	ZZ_pX		modulus;
	int 		size;
};




#endif
