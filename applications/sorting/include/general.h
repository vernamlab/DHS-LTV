#ifndef GENERAL_H_
#define GENERAL_H_

#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <time.h>

using namespace std;
NTL_CLIENT

//////////////////////////////////////////////////////////////////////////////////

typedef struct{
	ZZ q;
	ZZ B;
	int M_CycDegree;
	int N_PolyDegree;
	int D_FactDegree;
	int FactSize;
}GlobalParam;

void 	Set(GlobalParam &gp, ZZ q, ZZ B, int M);
ZZ 		euler_toient(ZZ x);
int		ComputeFactorDegree(int m);
void 	RandomPolyGen(ZZX &x, int N, int bs, ZZ q);

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

//////////////////////////////////////////////////////////////////////////////////

class myReduction{
public:
	void Set_degree(int n);
	void SetModulus(ZZX mod);
	void ComputeTable();
	void Reduce(ZZX &out, ZZX &in);
	void Div_high_low(ZZX &high, ZZX &low, ZZX &t, int low_deg, int high_deg);
private:
	ZZX modulus;
	int degree_n;
	ZZX *x_n;
	int loop;
};

//////////////////////////////////////////////////////////////////////////////////


#endif

















