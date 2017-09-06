#include "general.h"

void Set(GlobalParam &gp, ZZ q, ZZ B, int M){
	gp.q = q;
	gp.B = B;
	gp.M_CycDegree = M;

	gp.N_PolyDegree = to_long(euler_toient(to_ZZ(gp.M_CycDegree)));
	gp.D_FactDegree = ComputeFactorDegree(gp.M_CycDegree);
	gp.FactSize = gp.N_PolyDegree/gp.D_FactDegree;
}

ZZ euler_toient(ZZ x){
	if(x < to_ZZ("3"))
		return x;

	ZZ res = x;
	ZZ t = to_ZZ("2");
	bool is = false;

	while(x != to_ZZ("1")){
		while(GCD(x,t) == t){
			x=x/t;
			is = true;
		}
		if(is)
			res = res*(t-1)/t;
		is = false;
		t = NextPrime(t+1);
	}
	cout << "Phi of m is " << res << endl;
	return res;
}

int ComputeFactorDegree(int m){
	int p=1;
	bool loop = true;

	ZZ t;
	while(loop){
		t = power(to_ZZ("2"), p)-1;
		if(t%m == 0)
			loop = false;
		else
			p++;
	}
	return p;
}

void RandomPolyGen(ZZX &x, int N, int bs, ZZ q){
	clear(x);
	for(int i=0; i<N; i++)
		SetCoeff(x, i, RandomBits_ZZ(bs)%q);
}


//////////////////////////////////////////////////////////////////////////////////


void myTimer::Start(){
	startTime = clock();
}

void myTimer::Stop(){
	stopTime = clock();
}

double myTimer::GetTime(){
	return (double)(stopTime-startTime)/CLOCKS_PER_SEC;
}

void myTimer::ShowTime(string s){
	cout << s << ":\t" << GetTime() << endl;
}

///////////////////////////////////////////////////////////////////////////////////

void myReduction::SetModulus(ZZX mod){
	modulus = mod;
}

void myReduction::ComputeTable(){

	ZZX temp;
	for(int i=0; i<loop; i++){
		clear(temp);
		SetCoeff(temp, degree_n+(degree_n>>(i+1)), 1);
		x_n[i] = temp % modulus;
	}
}

void myReduction::Set_degree(int n){
	degree_n = n;
	loop = 0;
	while(n>128){
		loop++;
		n=n/2;
	}
	x_n = new ZZX[loop];
}

void myReduction::Div_high_low(ZZX &high, ZZX &low, ZZX &t, int low_deg, int high_deg){
	for(int i=0; i<low_deg; i++)
		SetCoeff(low, i, coeff(t, i));

	for(int i=low_deg; i<high_deg+1; i++)
		SetCoeff(high, i-low_deg, coeff(t, i));
}

void myReduction::Reduce(ZZX &out, ZZX &in){

	int low_deg, high_deg;
	ZZX high, low, t;
	t=in;
	for(int i=0; i<loop; i++){
		high_deg = degree_n+(degree_n>>(i));
		low_deg = degree_n+(degree_n>>(i+1));

		Div_high_low(high, low, t, low_deg, high_deg);
		t = low + high*x_n[i];
		clear(low);
		clear(high);
	}
	out = t%modulus;
}

///////////////////////////////////////////////////////////////////////////////////


















