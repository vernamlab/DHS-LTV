#include "ntt.h"

/*
BatcherNTT::BatcherNTT(GlobalParam *gp)
{
    Init(gp->p, gp->n);
    gp->p = to_int(p_);
}
*/
BatcherNTT::BatcherNTT(int p, int degree)
{
    Init(p, degree);
}


void BatcherNTT::Init(int p, int degree)
{
    n_ = degree;
    FindP(p);
    //p_ = p;
    //FindDegree();
    ZZ e;
	EulerToient(e, to_ZZ(p));
	FactorizeNumber(factors_, powers_, e);
	FindGenerator(g_, factors_, p_);
	FindNthRootOfUnity(alpha_, g_, p_, to_ZZ(n_));
	InvMod(inv_alpha_, alpha_, p_);
}

void BatcherNTT::FindP(int p)
{
	//ZZ t = to_ZZ(1)<<p_size;       // t = 2^k
	//p_ = n_*n_;
	//inv_n_ = InvMod(to_ZZ(n_), p_);

/*
	ZZ q;
	ZZ t = to_ZZ(p);
	t = NextPrime(t);
	t = t/n_;
	q = t*n_+1;

	while(ProbPrime(q, 100) != 1){
		t++;
		q = t*n_+1;
	}
	p_ = q;
	inv_n_ = InvMod(to_ZZ(n_), p_);
	//PrintP();
*/


/*

    p_ = n_ + 1;
    ZZ e;
	EulerToient(e, p_);
	cout << "phi(p) : " << e << endl;
*/

    // We want to find an n^th root of unity mod p eventually for FFT
    // Let it be w; s.t w^n = 1 mod p
    // We know that w^{phi(p)} = 1 mod p
    // But we want n to be the smallest such exponent
    // Then n | phi(p)
    // And n must be a power of two, for Cooley Tukey FFT
    // And 2 and p must be coprime, because we also
    // need to find (n-1)^{-1} mod p for FFT

    ZZ q;
    q = n_+1;

	while(ProbPrime(q, 100) != 1){
		q += n_;
	}
	p_ = q;

    //EulerToient(q, to_ZZ(n_));
    //n_ = to_int(q);
    //inv_n_ = InvMod(q, p_);
    inv_n_ = InvMod(to_ZZ(n_), p_);
    cout << "n : " << n_ << endl;
    cout << "n^{-1} : " << inv_n_ << endl;
}

void BatcherNTT::PrintFactors()
{
    ZZ e;
	EulerToient(e, p_);
    cout << e << " : ";
    for(int i=0; i<factors_.length(); i++)
    {
        if(i>0) cout << " + ";
        cout << factors_[i] << "^" << powers_[i];
    }
    cout << endl;
}

void BatcherNTT::fft_core(ZZX &out, const ZZX &in, const ZZ &r, int l)
{
	ZZX res;
	ZZ ta1, ta2;

	for(int i=0; i<l; i++)
	{
		ta1 = 0;
		for(int j=0; j<l; j++){
			ZZ tt = PowerMod(r, i*j,  p_);
			ta2 = coeff(in, j) * tt;
			ta1 = ta1 + ta2;
		}
		SetCoeff(res, i, ta1% p_);
	}

	out = res;
}

void BatcherNTT::CooleyTukey(ZZX &out, const ZZX &in, const ZZ &r, int l)
{
    ZZX X_FFT;
	ZZX X_even, X_odd;
	ZZX X_E_FFT, X_O_FFT;

	ZZ r2, t;

	X_FFT.SetLength(l);
	X_even.SetLength(l/2);
	X_odd.SetLength(l/2);
	X_E_FFT.SetLength(l/2);
	X_O_FFT.SetLength(l/2);

	if(l == 4)
        fft_core(X_FFT, in, r, l);

	else
	{
		for(int i=0; i<l/2; i++)
		{
			X_even[i] = coeff(in, 2*i);
			X_odd[i] = coeff(in, 2*i+1);
		}

		r2 = PowerMod(r, 2, p_);
		CooleyTukey(X_E_FFT, X_even, r2, l/2);
		CooleyTukey(X_O_FFT, X_odd,  r2, l/2);

		for(int i=0; i<l/2; i++)
		{
			t = PowerMod(r, i, p_);
			t = (t*X_O_FFT[i]) % p_;
			X_FFT[i] = (X_E_FFT[i] + t)% p_;
			X_FFT[i+l/2] = (X_E_FFT[i] - t)% p_;
		}
	}

	out = X_FFT;
}

void BatcherNTT::Encode(ZZX &out, const ZZX &in)
{
    CooleyTukey(out, in, inv_alpha_, n_);
	for(int i=0; i<=deg(out); i++)
		out[i] = (out[i]*inv_n_) % p_;

	//for(int i=0; i<n_; i++)
	//	out[i] = (out[i]*inv_n_) % p_;
}

void BatcherNTT::Decode(ZZX &out, const ZZX &in)
{
/*
	CooleyTukey(out, in, inv_alpha_, n_);
	for(int i=0; i<deg(out); i++)
		out[i] = (out[i]*inv_n_) % p_;
*/
    CooleyTukey(out, in, alpha_, n_);
}


////////////////////////////////////////////
////////////////////////////////////////////

/*
BatcherCRT::BatcherCRT(const GlobalParam *gp)
{
    m = gp->m;
    n = gp->n;
    d = gp->d;
    r = gp->r;
    p = gp->p;
    //////////////////////////////////////////////////
    mod_file = "CYC" + std::to_string(m) + ".txt";
    factor_file = "Factors_m" + std::to_string(m) + "_n" + std::to_string(n) + "_p" + std::to_string(p) +".txt";
    product_file = "CRT_MxN_m"+ std::to_string(m) +".txt";
}
*/
BatcherCRT::BatcherCRT(int cyc_m, int deg, int prime_p, int fac_deg, int fac_cnt)
{
    m = cyc_m;
    n = deg;
    d = fac_deg;
    r = fac_cnt;
    p = prime_p;
    //////////////////////////////////////////////////
    mod_file = "CYC" + std::to_string(m) + ".txt";
    factor_file = "Factors_m" + std::to_string(m) + "_n" + std::to_string(n) + "_p" + std::to_string(p) +".txt";
    product_file = "CRT_MxN_m"+ std::to_string(m) +".txt";
}
BatcherCRT::BatcherCRT(const ZZ &cyc_m, const ZZ &deg, const ZZ &prime_p, int fac_deg, int fac_cnt)
{
    m = to_int(cyc_m);
    n = to_int(deg);
    d = fac_deg;
    r = fac_cnt;
    p = to_int(prime_p);
    //////////////////////////////////////////////////
    mod_file = "CYC" + std::to_string(m) + ".txt";
    factor_file = "Factors_m" + std::to_string(m) + "_n" + std::to_string(n) + "_p" + std::to_string(p) +".txt";
    product_file = "CRT_MxN_m"+ std::to_string(m) +".txt";
}

void BatcherCRT::Encode(ZZX &out, const ZZX &in)
{
    ZZ_p::init(to_ZZ(p));
    ifstream file;
    file.open(product_file);
	ZZ_pX res, temp, mod;
	res = 0;
	for(int i=0; i<r && i<=deg(in); i++)
	{
        file >> temp;
        temp = to_ZZ_pX(coeff(in,i)*to_ZZX(temp));
        res += temp;
    }
    file.close();
    file.open(mod_file);
    file >> mod;
	res = res % mod;
	out = to_ZZX(res);
	file.close();
}

void BatcherCRT::Decode(ZZX &out, const ZZX &in)
{
    ZZ_p::init(to_ZZ(p));
    ifstream file;
    file.open(factor_file);
    ZZ_pX factor, msg, temp;
    msg = to_ZZ_pX(in);
    for(int i=0; i<r && i<deg(in); i++)
    {
        file >> factor;
        temp = msg%factor;
        ZZ t = rep(coeff(temp,0));
        SetCoeff(out, i, t%to_ZZ(p));
    }
    file.close();
}
