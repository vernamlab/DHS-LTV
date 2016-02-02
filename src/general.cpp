#include "general.h"
/*
typedef struct{
    // Application Params
    int p;      // plaintext space, 2 or 257?
    int depth;

    // Security Params
    int k;      // cutting factor
    int l;      // how many times we will cut?
    int i;      // the last level q size
	ZZ* q_list; // q[i] = k^i
	ZZ B;       // noise bound, usually 1

	// Ring Params
	int m;  // cyc m
	int n;  // cyc degree n = Phi(m)
	int d;  // fact degree
	int r;  // fact cnt n = d.r, slot count
}GlobalParam;
*/

// Operator Overloads
CipherText& CipherText::operator=(const CipherText& rhs)
{
    if(this != &rhs) // If its not a self assignment
    {
        ct_ = rhs.ct();
        level_ = rhs.level();
        fresh_ = rhs.fresh();
        multcnt_ = rhs.multcnt();
    }
    return *this;
}
CipherText& CipherText::operator+=(const CipherText& rhs)
{
    set_ct(ct()+rhs.ct());
    return *this;
}
const CipherText CipherText::operator+(const CipherText& rhs) const
{
    CipherText ct;
    ct = *this;
    ct += rhs;
    return ct;
}
CipherText& CipherText::operator-=(const CipherText& rhs)
{
    set_ct(ct()-rhs.ct());
    return *this;
}
const CipherText CipherText::operator-(const CipherText& rhs) const
{
    CipherText ct;
    ct = *this;
    ct -= rhs;
    return ct;
}
CipherText& CipherText::operator+=(const ZZ& rhs)
{
    set_ct(ct()+rhs);
    return *this;
}
const CipherText CipherText::operator+(const ZZ& rhs) const
{
    CipherText ct;
    ct = *this;
    ct += rhs;
    return ct;
}
CipherText& CipherText::operator*=(const ZZ& rhs)
{
    set_ct(ct()*rhs);
    return *this;
}
const CipherText CipherText::operator*(const ZZ& rhs) const
{
    CipherText ct;
    ct = *this;
    ct *= rhs;
    return ct;
}
std::ostream& operator<<(std::ostream& os, const CipherText& ct)
{
    os  << ct.ct();
    return os;
}

void CipherText::coeff_reduce(const ZZ &mod)
{
    CoeffReduce(ct_, ct_, mod);
}
void CipherText::poly_reduce(RingType type, const ZZX &mod)
{
    PolyReduce(ct_, ct_, type, mod);
}

GlobalParams::GlobalParams()
{
    set_l(1);   // Circuit Level
    set_hom_type(she);
    set_reduc_type(yarkin);
    Init(cyclotomic, off, 2, 8191);
}

GlobalParams::GlobalParams(RingType type, BatchFlag flag, int p, int degree)
{
    Init(type, flag, p, degree);
}
void GlobalParams::Init(RingType type, BatchFlag flag, int p, int degree)
{
    set_ring_type(type);
    set_batch_flag(flag);
    if(type == xn_1)
    {
        set_n(degree);
        if(flag == on)
        {
            set_d(1);
            set_p(FindPrimeCongOne(n_));    // p = 1 mod n
        }
        else
            set_p(p);
    }
    else if(type == cyclotomic)
    {
        if(flag == on)
        {
            if(p == 2)
            {
                set_m(degree);
                set_d(ComputeFactorDegree(m_));
            }
            else
            {
                set_d(2);
                set_m(p*p-1);   // m = p^2 - 1
            }
        }
        else
            set_m(degree);

        set_n(EulerToient(m_));
        set_p(p);
    }

    // Security Params
    set_b(1);
    set_k(100); // (w*2)*(l+1)

    if(hom_type_ == fhe)
        set_i(k_);
    else
        set_i(100);

    //set_l(1);

    if(flag == on)
        set_r(n_/d_);
    else
        set_r(1);

    set_w(20);   // How many bits of secret key will be encrypted
}

/*
void FindQList(vec_ZZ &q_list, const GlobalParam *gp)
{
    q_list.SetLength(gp->l+1);
	ZZ 	t, co;
	GenPrime(t, gp->i);
	while(t%gp->p != 1)
        GenPrime(t, gp->i);

	q_list[gp->l] = t;//Smallest q, last one

    GenPrime(t, gp->k);
    while(t%gp->p != 1)
        GenPrime(t, gp->k);

	for(int i=gp->l-1; i>=0; i--)
		q_list[i] = q_list[i+1]*t;
}
*/


/****************************************************************************************/
/************************************ POLY OPERATIONS ***********************************/
/****************************************************************************************/
void RandPolyBalanced(ZZX &out, int n, int q)
{
    RandPolyBalanced(out, n, to_ZZ(q));
}

void RandPolyBalanced(ZZX &out, int n, const ZZ &q)
{
    // deg(out) = n-1
    // -q <= out[i] <= q
    //clear(out);
    //cout << n << endl;
	out.SetLength(n);
	for(int i=0; i<n; i++)
	{
		out[i] = RandomBnd(2*q+1);
        out[i] -= q;
    }
    out.normalize();
}

void RandPolyBounded(ZZX &out, int n, int q)
{
    RandPolyBounded(out, n, to_ZZ(q));
}

void RandPolyBounded(ZZX &out, int n, const ZZ &q)
{
    // deg(out) = n-1
    // 0 <= out[i] <= q-1
    //clear(out);
	out.SetLength(n);
	for(int i=0; i<n; i++)
		out[i] = RandomBnd(q);
    out.normalize();
}

void CyclotomicPoly(ZZX &out, int m)
{
    //clear(out);
	int s;
	out = 1;
	int n = m;

	ZZX t_vec[n];
	int s_vec[n];

	for(int i=0; i<n; i++)
		s_vec[i] = 0;

	for(int d=1; d<=n; d++){
			if(GCD(d, n) == d){

				ZZX t;
				SetCoeff(t, 0 , -1);
				SetCoeff(t, n/d, 1);

				s = MobuisFunction(d);

				t_vec[d-1] = t;
				s_vec[d-1] = s;
			}
	}

	for(int i=0; i<n; i++)
		if(s_vec[i] == 1)
			out = out * t_vec[i];

	for(int i=0; i<n; i++)
		if(s_vec[i] == -1)
			out = out /  t_vec[i];
}

void XN_1(ZZX &out, int n)
{
    SetCoeff(out, n, 1);
    SetCoeff(out, 0, -1);
}

/**
    Polynomial U;
    U = x^{2*n-1};
    U /= mod;
**/
void BarrettReduction(ZZX &out, const ZZX &in, const ZZX &mod, const ZZX &u)
{
    ZZX t = in;
    int n = deg(mod);
    ZZX m = mod;
    ZZX X_N;
    SetCoeff(X_N, n, 1);
    ////////////////////
    ZZX q;
    div(q, t, X_N);             // 2
    q *= u;                     // 3
    X_N = X_N >> 1;
    div(q, q, X_N);             // 4
    q *= m;                     // 5
    q = t - q;                  // 6
    while(deg(q) >= deg(m))        // 7
        q -= m;
    out = q;
}
void BarrettReduction(ZZX &out, const ZZX &in, const ZZX &mod)
{
    ZZX t = in;
    int n = deg(mod);

    ZZX u, r;
    SetCoeff(u, 2*n-1, 1);
    u /= mod;

    ZZX m = mod;

    ZZX X_N;
    SetCoeff(X_N, n, 1);

    //////////////////////////////
    ZZX q;
    div(q, t, X_N);             // 2
    q *= u;                     // 3
    X_N = X_N >> 1;
    div(q, q, X_N);             // 4
    q *= m;                     // 5
    q = t - q;                  // 6
    while(deg(q) >= deg(m))        // 7
        q -= m;
    out = q;
}
void Div_high_low(ZZX &high, ZZX &low, ZZX &t, int low_deg, int high_deg)
{
	for(int i=0; i<low_deg; i++)
		SetCoeff(low, i, coeff(t, i));

	for(int i=low_deg; i<high_deg+1; i++)
		SetCoeff(high, i-low_deg, coeff(t, i));
}
void  YarkinReduction(ZZX &out, const ZZX &in, const ZZX &mod, const vec_ZZX &helpers)
{
    int low_deg, high_deg;
	ZZX high, low, t;
	t = in;
	int n = deg(mod);
	int cnt = helpers.length();
	for(int i=0; i<cnt; i++)
	{
		high_deg = n + (n>>i);
		low_deg = n + (n>>(i+1));

		Div_high_low(high, low, t, low_deg, high_deg);
		t = low + high * helpers[i];
		clear(low);
		clear(high);
	}
	out = t % mod;
}
void YarkinReduction(ZZX &out, const ZZX &in, const ZZX &mod)
{
    int cnt = 0;
    int n = deg(mod);
    int temp = n;
	while(temp > 128)
	{
		cnt++;
		temp /= 2;
	}

	vec_ZZX polylist;
	polylist.SetLength(cnt);

	for(int i=0; i<cnt; i++)
	{
		SetCoeff(polylist[i], n+(n>>(i+1)), 1);
		polylist[i] = polylist[i] % mod;
	}

	int low_deg, high_deg;
	ZZX high, low, t;
	t = in;
	for(int i=0; i<cnt; i++)
	{
		high_deg = n + (n>>i);
		low_deg = n + (n>>(i+1));

		Div_high_low(high, low, t, low_deg, high_deg);
		t = low + high * polylist[i];
		clear(low);
		clear(high);
	}
	out = t % mod;
}

void PolyReduce(ZZX &out, const ZZX &in, RingType type, const ZZX &mod)
{
    out = in % mod;
	out.normalize();
}

void CoeffReduce(ZZX &out, const ZZX &in, int q)
{
    CoeffReduce(out, in, to_ZZ(q));
}

void CoeffReduce(ZZX &out, const ZZX &in, const ZZ &mod)
{
    out = in;
    for(int i=0; i<=deg(in); i++)
        out[i] = out[i] % mod;
    out.normalize();
}

void CoeffBalance(ZZX &out, const ZZX &in, const ZZ &q)
{
    out = in;
    for(int i=0; i<=deg(in); i++)
    {
        out[i] = out[i] % q;
        if(out[i] > ((q-1)/2))
			out[i] -= q;
    }
    out.normalize();
}

void CoeffDivide(ZZX &out, const ZZX &in, const ZZ &q)
{
    out = in;
    for(int i=0; i<=deg(in); i++)
    {
        out[i] = out[i] / q;
    }
    out.normalize();
}

void CoeffParityMatch(ZZX &out, const ZZX &in, int q)
{
    out = in;
    for(int i=0; i<=deg(in); i++)
    {
        out[i] = out[i] - (out[i]%q) + (in[i]%q);
    }
    out.normalize();
}

void PolyCoeffReduce(ZZX &out, const ZZX &in, RingType type, const ZZX &poly_mod, const ZZ &coeff_mod)
{
    PolyReduce(out, in, type, poly_mod);
    CoeffReduce(out, out, coeff_mod);
}
void Clear(ZZX &in)
{
	for(int i=0; i<deg(in); i++)
		in[i] = 0;
}
void FindInverse(ZZX &out, const ZZX &in, const ZZX& poly_mod, const ZZ &coeff_mod, bool &found)
{
	ZZ_p::init(coeff_mod);
	ZZ_pX phi;
	phi = to_ZZ_pX(poly_mod);
	ZZ_pE::init(phi);

	ZZ_pE f_, f_inv_;
	f_ = to_ZZ_pE(to_ZZ_pX(in));
	try{
        f_inv_ = inv(f_);
        ZZ_pX tp = rep(f_inv_);
        out = to_ZZX(tp);
        found = true;
	}
	catch(int e){
        cout << "Error : could not find an inverse." << endl;
        found = false;
    }
}
bool FindInverse(ZZX &out, const ZZX &in, const ZZX& poly_mod, const ZZ &coeff_mod)
{
	bool result;
	FindInverse(out, in, poly_mod, coeff_mod, result);
	return result;
}

void ComputeFFT(fftRep *rep_list, const ZZX &poly, int prime_cnt, int deg)
{
    zz_pBak bak;
    bak.save();
    zz_pX temp;
    long k = NextPowerOfTwo(2*deg+1);
    for(int i=0; i<prime_cnt; i++)
    {
        zz_p::FFTInit(i);
        conv(temp, poly);
        TofftRep(rep_list[i], temp, k);
    }
    bak.restore();
}

void FFTRepClear(fftRep *rep_list, int length)
{
    long* p;
	for(int i=0; i<length; i++)
	{
		p = &(rep_list[i].tbl[0][0]);
		for(int j=0; j<(1L<<rep_list[i].k); j++)
			p[j] = 0;
	}
}

void FFTRepSetK(fftRep *rep_list, int length, int k)
{
    for(int i=0; i<length; i++)
    {
		zz_p::FFTInit(i);
		rep_list[i].SetSize(k);
	}
}

void SeparateBlocks(vec_ZZX &out, const ZZX &in, long blocksize)
{
    ZZ temp = to_ZZ(1);
    temp = temp<<blocksize;    // 2^w
    ZZX msg = in;
    for(int i=0; i<out.length(); i++)
    {
        out[i].SetLength(deg(msg)+1);
		for(int j=0; j<=deg(msg); j++)
        {
            out[i][j] = msg[j]%temp;
            msg[j] = msg[j]/temp;
        }
    }
}

void FFTtoPoly(ZZX &out, vec_zz_pX &in, int rep_length, int n)
{

    zz_pBak bak;
    bak.save();

    ZZ prod;
    set(prod);
	for(int i=0; i<rep_length; i++)
		mul(prod, prod, GetFFTPrime(i));

    ZZ temp1, temp2, res;
    vec_ZZ c;
    c.SetLength(2*n+1);
    for(int i = 0; i<rep_length; i++)
    {
		zz_p::FFTInit(i);
	    long p = zz_p::modulus();

	    temp1 = prod / p;
	    temp2 = temp1 % p;
	    temp2 = InvMod(temp2, to_ZZ(p));
	    res = temp1 * temp2;

	    long m = deg(in[i]);

	    for(int j=0; j<=m; j++)
	    {
	    	mul(temp1, res, rep(in[i].rep[j]));
	        add(c[j], c[j], temp1);
        }
	}


    out.rep.SetLength(2*n+1);
	ZZ prod2 = prod >> 1;

	//RightShift(prod2, prod, 1);

	for(int j=0; j<=2*n; j++)
	{
        temp1 = c[j] % prod;
	    if (temp1 > prod2)
	    	out.rep[j] = temp1 - prod;
	    else
	    	out.rep[j] = temp1;
	}
	   out.normalize();
	   bak.restore();

/*
void ntru::FromCrt2Poly(ZZX &res, zz_pX *polys, int &priNum){
	 zz_pBak bak;
	 bak.save();

	ZZ prod, coeff, t1;
	long tt;
	vec_ZZ c;

	c.SetLength(N+N+1);

	set(prod);
	for(int i=0; i<priNum; i++)
		mul(prod, prod, GetFFTPrime(i));

	for(int i = 0; i < priNum; i++) {
		zz_p::FFTInit(i);
	    long p = zz_p::modulus();

	    div(t1, prod, p);
	    tt = rem(t1, p);
	    tt = InvMod(tt, p);
	    mul(coeff, t1, tt);

	    long m = deg(polys[i]);

	    for(int j = 0; j <= m; j++){
	    	mul(t1, coeff, rep(polys[i].rep[j]));
	        add(c[j], c[j], t1);
	      }
	}

	res.rep.SetLength(N+N+1);

	ZZ prod2;
	RightShift(prod2, prod, 1);

	for(int j = 0; j <= N+N; j++){
		rem(t1, c[j], prod);
	    if (t1 > prod2)
	    	sub(res.rep[j], t1, prod);
	    else
	    	res.rep[j] = t1;
	}
	   res.normalize();
	   bak.restore();
}
*/
}


/******************************************************************************************/
/********************************** BIG NUMBER OPERATIONS *********************************/
/******************************************************************************************/

long GetPrimeNumber(long bound, ZZ &prod){

	long nprimes;
	zz_pBak bak;
	bak.save();
	//prod = 1;
	for (nprimes = 0; NumBits(prod) <= bound; nprimes++) {
			UseFFTPrime(nprimes);
		mul(prod, prod, GetFFTPrime(nprimes));
   }
	bak.restore();
	return nprimes;
}

void FindPrimeCongOne(ZZ &out, const ZZ &mod)
{
    out = mod+1;

	while(ProbPrime(out, 100) != 1){
		out += mod;
	}
}

int FindPrimeCongOne(int mod)
{
    ZZ result;
    FindPrimeCongOne(result, to_ZZ(mod));
    return to_int(result);
}

void EulerToient(ZZ &out, const ZZ &in)
{
    out = in;
	if(in > to_ZZ("2"))
    {
        ZZ x = in;
        ZZ t = to_ZZ("2");
        bool is = false;

        while(x != to_ZZ("1")){
            while(GCD(x,t) == t){
                x=x/t;
                is = true;
            }
            if(is)
                out = out*(t-1)/t;
            is = false;
            t = NextPrime(t+1);
        }
	}
}
int EulerToient(int in)
{
    ZZ result;
    EulerToient(result, to_ZZ(in));
    return to_int(result);
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
/**
*** f = f_1^k_1 . f_2^k_2 . . . f_t^k_t
**/
void FactorizeNumber(vec_ZZ &factors, vec_ZZ &powers, const ZZ &in)
{
    ZZ f = in;
    ZZ s = to_ZZ("2");
	int l = 0;
	while(s<=f && f>1)
	{
		if(f%s == 0)
		{
            factors.SetLength(l+1);
            powers.SetLength(l+1);

			factors[l] = s;
			powers[l] = 1;

			f = f/s;
			while(f%s == 0)
			{
                f = f/s;
                powers[l]++;
			}
			l++;
		}
		s++;
		s = NextPrime(s);
	}
}
void FindGenerator(ZZ &out, const vec_ZZ &factors, const ZZ &p)
{
    ZZ f;
    EulerToient(f, p);
	ZZ s = to_ZZ(1);

	bool found = false;
	while(!found)
	{
        s++;
        found = true;
        for(int i=0; i<factors.length() && found; i++)
        {
            ZZ exp = f/factors[i];
            ZZ temp = PowerMod(s, exp, p);
            if(temp == 1)
                found = false;
        }
    }
	out = s;
}
void FindNthRootOfUnity(ZZ &out, const ZZ &g, const ZZ &p, const ZZ &n)
{
    ZZ f;
    EulerToient(f, p);
    ZZ t = f/n;
	PowerMod(out, g, t, p);
}
int MobuisFunction(int n){
	int t, primes;
	primes = 0;

	if(n == 1)
		return 1;
	else{
		for(int i=2; i<=n; i++){
			if(ProbPrime(i)){
				if(GCD(i,n) == i){
					t=n/i;
					primes++;
					if(GCD(i, t) == i)
						return 0;
				}
			}
		}
		if(primes%2 == 0)
			return 1;
		else
			return -1;
	}
}

/*************************************************************************************/
/**************************************** TIMER **************************************/
/*************************************************************************************/

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


















