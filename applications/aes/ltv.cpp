#include "ltv.h"

//////////////////////////////////////////////////////////////////////////////////////////////
ltv::ltv(int dm, int num, int max_bit, int tableS, int wordSize){
	word 	= wordSize;
	ZZ qq 	= to_ZZ("2");
	ZZ BB 	= to_ZZ("1");
	Set(gp, qq, BB, dm);

	p 			= gp.q;
	B 			= gp.B;
	N 			= gp.N_PolyDegree;
	cout << "N_PolyDegree " << N << endl;
	degree_m 	= gp.M_CycDegree;
	cout << "CycDegree " << degree_m << endl;
	rand_seed 	 = 0;
	message_rand = 0;
	reducsetflag = false;

	clear(modulus, N+1);
	ZZ_p::init(p);


	tableSize = tableS;
	word_blocksize = NumofWords(max_bit);
	max_bitsize = max_bit;
	tableSize 	 = tableS;

	pk = new ZZX[num];
	sk  = new ZZX[num];
	q_list = new ZZ[num];
	ek2 = new struct eval_key;
	ek2[0].key 	= new ZZX[word_blocksize];


	fftKeys = new FFTPolyList[tableSize];
	for(int i=0; i<tableSize; i++){
		fftKeys[i].polynum		= word_blocksize;
		fftKeys[i].polyList 	= new PolyFFTRep[word_blocksize];
		for(int j=0; j<word_blocksize; j++)
			fftKeys[i].polyList[j].poly_rep = new fftRep[5];
	}


	//Set Randomness
#ifdef PureRandom
	srand(time(NULL));
	SetSeed(to_ZZ(time(NULL)));
#endif

#ifdef PseudoRandom
	srand(5);
	SetSeed(to_ZZ("5"));
#endif
}

ltv::ltv(ZZ my_p, ZZ my_B, int my_N, int dm, int num, int max_bit, int tableS, int wordSize){
	word 	= wordSize;
	p = my_p;	B = my_B;
	N = my_N;	degree_m = dm;
	rand_seed 	 = 0;
	message_rand = 0;
	reducsetflag = false;

	clear(modulus, N+1);
	ZZ_p::init(p);

	tableSize = tableS;
	word_blocksize = NumofWords(max_bit);	// maxprime/wordsize+1
	max_bitsize = max_bit;
	tableSize 	 = tableS;

	pk = new ZZX[num];
	sk  = new ZZX[num];
	q_list = new ZZ[num];
	ek2 = new struct eval_key;
	ek2[0].key 	= new ZZX[word_blocksize];

	fftKeys = new FFTPolyList[tableSize];
	for(int i=0; i<tableSize; i++){
		fftKeys[i].polynum		= word_blocksize;
		fftKeys[i].polyList 	= new PolyFFTRep[word_blocksize];
		for(int j=0; j<word_blocksize; j++)
			fftKeys[i].polyList[j].poly_rep = new fftRep[5];
	}
}

int ltv::NumofWords(int bs){
	int size;

	size = bs/word;
	if (bs%word != 0)
		size = size+1;

	return size;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


ZZX ltv::Prim_Encrypt(ZZX m, int i){

	ZZX s, e, result;
	s = Func_Sample();
	e = Func_Sample();
	coeff_reduction(s, s, q_list[i], N);
	coeff_reduction(e, e, q_list[i], N);
	result = pk[i]*s + e*p + m;

	Arith_PolyReduce(result, result);
	coeff_reduction(result, result, q_list[i], N);

	return result;
}

ZZX ltv::Prim_Decrypt(ZZX c, int index){

	ZZ x;
	ZZX result, t;

	t = sk[index]*c;
	Arith_PolyReduce(t,t);
	coeff_reduction(t, t, q_list[index], N);
	for(int i=0; i<N; i++){
		x = coeff(t, i);
		if(x>((q_list[index]-1)/2))
			x = x-q_list[index];
		SetCoeff(result, i, (x%p));
	}
	return result;
}


ZZX ltv::Prim_RelinRingFFT(ZZX &c, int index, int tableIndx){

	zz_pBak bak;
	bak.save();
	int priNum = fftKeys[tableIndx].polyList[0].rep_num;

	ZZX bits, res;
	fftRep *fft_bit, *fft_res, *fft_mul;
	long fftSize = NextPowerOfTwo(2*N-1);

	fft_bit = new fftRep [priNum];
	fft_res = new fftRep [priNum];
	fft_mul = new fftRep [priNum];

	fftRepSetSize(fft_bit, priNum, fftSize);
	fftRepSetSize(fft_res, priNum, fftSize);
	fftRepSetSize(fft_mul, priNum, fftSize);

	fftRepClear(fft_bit, priNum);
	fftRepClear(fft_res, priNum);
	fftRepClear(fft_mul, priNum);

	ZZX rt;
	int leftWords = max_bitsize-index*Dif_Prime;
	leftWords = NumofWords(leftWords);

	for(int j=0; j<(leftWords); j++){
//		GetPolyBitsIndex(bits, c, j);
		GetPolyWordsIndex(bits, c, j);
		CalculateFFTValues(fft_bit, bits, priNum, N-1);
		FFTmult(fft_mul, fft_bit, fftKeys[tableIndx].polyList[j].poly_rep, priNum);
		FFTadd(fft_res, fft_res, fft_mul, priNum);
	}

	zz_pX rescrt[priNum];
	for(int i=0; i<priNum; i++){
		zz_p::FFTInit(i);
		FromfftRep(rescrt[i], fft_res[i], 0, N+N-2);
	}
	FromCrt2Poly(res, rescrt, priNum);

	int y = 2*N;
	coeff_reduction(res, res, q_list[index], y);

	delete [] fft_bit;
	delete [] fft_res;
	delete [] fft_mul;

	bak.restore();
	return res;
}

ZZX ltv::Prim_RelinRingFFT(ZZX &c, int index){
	return Prim_RelinRingFFT(c, index, 0);
}


ZZX ltv::Prim_ModSwitch(ZZX &x, int index){
	ZZX res;
	ZZ tp;
	for(int i=0; i<N; i++){
		tp = coeff(x, i);
		if(tp>((q_list[index]-1)/2))
			tp = tp - q_list[index];

		getClosestMod2(tp, tp, q_list[index], q_list[index+1]);
		SetCoeff(res, i, tp);
	}
	return res;
}

ZZX ltv::Prim_ModSwitch_LongPoly(ZZX &x, int PolyDeg, int index){
	ZZX res;
	ZZ tp;
	for(int i=0; i<PolyDeg; i++){
		tp = coeff(x, i);
		if(tp>((q_list[index]-1)/2))
			tp = tp - q_list[index];

		getClosestMod2(tp, tp, q_list[index], q_list[index+1]);
		SetCoeff(res, i, tp);
	}
	return res;
}

void ltv::Prim_ReduceKeysLevel(int level){
	Prim_ReduceKeysLevel(level, 0);
}

void ltv::Prim_ReduceKeysLevel(int level, int tblIndex){

	int leftBits = max_bitsize-level*Dif_Prime;
	int leftBlocks = NumofWords(leftBits);

	for(int i=0; i<leftBlocks; i++)
		coeff_reduction(ek2[0].key[i], ek2[0].key[i], q_list[level], N);

//	for(int i=0; i<leftBits; i++){
//		ek2[0].key[i] = Prim_ModSwitch(ek2[0].key[i], level);
//		coeff_reduction(ek2[0].key[i], ek2[0].key[i], q_list[level], N);
//	}

	for(int i=leftBlocks; i<word_blocksize; i++)
		ek2[0].key[i] = 0;


	for(int i=0; i<fftKeys[tblIndex].polynum; i++)
		delete [] fftKeys[tblIndex].polyList[i].poly_rep;
	delete [] fftKeys[tblIndex].polyList;

	fftKeys[tblIndex].polynum		= leftBlocks;
	fftKeys[tblIndex].polyList 	= new PolyFFTRep[leftBlocks];

	ConvertKeystoFFT(ek2[0], fftKeys[tblIndex], max_bitsize, leftBlocks);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

ZZX ltv::Func_Sample(){
	ZZX result;
	for(int i = 0; i < N;i ++)
	{
		ZZ tp;
		RandomBits(tp,to_long(B));
		//cout << "RandomBits " << tp << " ";
		if(rand()%2)SetCoeff(result,i, tp);
		else SetCoeff(result,i, 0-tp);
	}
	//cout << endl;
	//cout << "random result " << result << endl;
	return result;
}

void ltv::Func_ModulusFindRing(int num, int max_bit, int diff, ZZX modu){

	ZZ 	t, co;
	max_bitsize = max_bit;

	GenPrime(t, diff);
	co = t;
	cout << "prime " << t << endl;
	for(int i=0; i<num; i++){
		q_list[num-1-i] = t;
		t=t*co;
	}

/*
	ZZ 	t, co;
	max_bitsize = max_bit;
	init_bitsize = max_bitsize;

	GenPrime(t, diff);
	co = t;
	for(int i=0; i<num; i++){
		q_list[num-1-i] = t;
		GenPrime(co, diff);
		cout << co << endl;
		t=t*co;
	}
*/

}

ZZX ltv::Func_CreateMessage(int size){
	ZZX x;
	ZZ r = to_ZZ(time(NULL));
	message_rand++;
	SetSeed(r+message_rand);

	long l;
	for(int i=0; i<size; i++){
		l = RandomBnd(2);
		SetCoeff(x, i, l);
	}
	return x;
}

void ltv::Func_SetModulus(){
	modulus = ComputeFastCycModulus(degree_m);
	cout << "modulus " << modulus << endl;
}

void ltv::Func_ComputeKeysRingRelinFFT(int num, int tblsize){
	tableSize = tblsize;

	if(reducsetflag == false)
		SetFastPolyModulus();

////////////////////////First Key//////////////////////////////////////
	ZZX f, g, ft, f_inv;
	bool isfound = false;

	while(!isfound){
		isfound = true;
		ft  = Func_Sample();
		coeff_reduction(ft, ft, q_list[0], N);
		f = ft*p + 1;
		coeff_reduction(f, f, q_list[0], N);
		find_inverse(f_inv, f, q_list[0], N, isfound, modulus);
		coeff_reduction(f_inv, f_inv, q_list[0], N);
	}
	isfound = false;
	g 	  = Func_Sample();
	coeff_reduction(g, g, q_list[0], N);
	pk[0] = g*f_inv;
	Arith_PolyReduce(pk[0], pk[0]);

	pk[0] = pk[0]*p;
	sk[0] = f;

	coeff_reduction(pk[0], pk[0], q_list[0], N);
	coeff_reduction(sk[0], sk[0], q_list[0], N);

	ZZX	tppk = pk[0];
	ZZX tpek = sk[0];

	compute_eval_onerelin(ek2[0], tppk, tpek, 0);
/////////////////////////

	ConvertKeystoFFT(ek2[0], fftKeys[0], max_bitsize, word_blocksize);

	for(int i=1; i<tblsize; i++)
		Prim_ReduceKeysLevel(i, i);

//////////////////////////Rest Key///////////////////////////////////
	for(int i=1; i<num; i++){
		sk[i] = sk[i-1];
		coeff_reduction(sk[i], sk[i], q_list[i], N);
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


void ltv::Arith_PolyReduce(ZZX &out, ZZX &in){
	myr.Reduce(out, in);
}

void ltv::Arith_CoeffReduce(ZZX &out, ZZX &in, int index){
	int l = in.rep.length();
	coeff_reduction(out, in, q_list[index], l);
}

void ltv::Arith_PolyCoeffReduce(ZZX &out, ZZX &in, int index){
	Arith_PolyReduce(out, in);
	Arith_CoeffReduce(out, out, index);
}

ZZX ltv::Arith_MulModZZX(ZZX &a, ZZX &b, int index){
	ZZX t;
	t = a*b;
	Arith_PolyReduce(t, t);
	coeff_reduction(t, t, q_list[index], N);
	return t;
}

ZZX ltv::Arith_AddModZZX(ZZX &a, ZZX &b, int index){
	ZZX t = (a+b);
	Arith_PolyReduce(t, t);
	coeff_reduction(t, t, q_list[index], N);
	return t;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

ZZ *ltv::Pointer_Q(int i){
	return &(q_list[i]);
}

ZZX	*ltv::Pointer_PolyModulus(){
	return (&modulus);
}

GlobalParam *ltv::Pointer_GP(){
	return (&gp);
}

ZZX	*ltv::Pointer_PublicKey(){
	return &(pk[0]);
}

ZZX	*ltv::Pointer_SecretKey(int i){
	return &(sk[i]);
}

ZZX	*ltv::Pointer_EvalKey(int i){
	return &(ek2[0].key[i]);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

ZZX ltv::ComputeFastCycModulus(int n){
	ZZX modulus;
	int s;
	modulus = 1;

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
			modulus = modulus * t_vec[i];

	for(int i=0; i<n; i++)
		if(s_vec[i] == -1)
			modulus = modulus /  t_vec[i];

	return modulus;
}

void ltv::SetFastPolyModulus(){
	myr.SetModulus(modulus);
	myr.Set_degree(N);
	myr.ComputeTable();
	reducsetflag = true;
}


void ltv::compute_eval_onerelin(eval_key &ek2, ZZX tppk, ZZX tpek, int index){

	ZZX tpek2 = tpek;

	ZZ tw =to_ZZ("1");
	ZZ w = to_ZZ("1");
	w = w << word;
	ZZX s, e, result, tp;
	for(int i=0; i<word_blocksize; i++){

		tp = tpek2 * tw;
		s = Func_Sample();
		e = Func_Sample();
		coeff_reduction(s, s, q_list[index], N);
		coeff_reduction(e, e, q_list[index], N);
		coeff_reduction(tp, tp, q_list[index], N);
		ek2.key[i] = tppk*s + e*p + tp;
		Arith_PolyReduce(ek2.key[i], ek2.key[i]);
		coeff_reduction(ek2.key[i], ek2.key[i], q_list[index], N);
		tw = tw*w;
	}
}


void ltv::ConvertKeystoFFT(eval_key &ek2, FFTPolyList &fftkey, int &bitsize, int &num_add){

	long bound 		= GetPrimeBound(N, N, bitsize, word, num_add);

	ZZ prod;
	set(prod);
	long priNum 	= GetPrimeNumber(bound, prod);
	zz_p::FFTInit(0);
	for(int i=0; i<fftkey.polynum; i++){
		fftkey.polyList[i].poly_rep = new fftRep[priNum];
		fftkey.polyList[i].rep_num	= priNum;
	}

	long mk = NextPowerOfTwo(2*N-1);
	for(int i=0; i<fftkey.polynum; i++)
		for(int j=0; j<fftkey.polyList[i].rep_num; j++){
			zz_p::FFTInit(j);
			fftkey.polyList[i].poly_rep[j].SetSize(mk);
		}

	for(int i=0; i<fftkey.polynum; i++)
		CalculateFFTValues(fftkey.polyList[i].poly_rep, ek2.key[i], fftkey.polyList[i].rep_num, N-1);

}

void ltv::GetPolyWordsIndex(ZZX &bits, ZZX &cip, int &bitIndex){

	ZZX t[word];
	int begin = bitIndex*word;

	for(int i=0; i<word; i++){
		GetPolyBitsIndex(t[i], cip, begin);

		for(int j=0; j<N; j++)
			SetCoeff(t[i], j, coeff(t[i], j)<<i);

		begin++;
	}

	ZZX r = t[0];
	for(int i=1; i<word; i++){
		r = r + t[i];
	}
	bits = r;
}


void ltv::GetPolyBitsIndex(ZZX &bits, ZZX &cip, int &bitIndex){
	for(int i=0; i<N; i++)
		SetCoeff(bits, i, bit(coeff(cip,i),bitIndex));
}

void ltv::FFTmult(fftRep *res, fftRep *R1, fftRep *R2, int &priNum){
	for(int i=0; i<priNum; i++){
		zz_p::FFTInit(i);
		mul(res[i], R1[i], R2[i]);
	}
}

void ltv::FFTadd(fftRep *res, fftRep *R1, fftRep *R2, int &priNum){
	for(int i=0; i<priNum; i++){
		zz_p::FFTInit(i);
		add(res[i], R1[i], R2[i]);
	}
}

void ltv::FromCrt2Poly(ZZX &res, zz_pX *polys, int &priNum){
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void clear(ZZX &x, int size){
	for(int i=0; i<size; i++)
		SetCoeff(x, i, to_ZZ("0"));
}

void coeff_reduction(ZZX &out, ZZX &in, ZZ &q, int &degree){
	for(int i=0; i<degree; i++)
		SetCoeff(out, i, (coeff(in,i)%q));
}

void find_inverse(ZZX &f_inv, ZZX &f, ZZ &q, int &degree, bool &isfound, ZZX modulus){

	ZZ_p::init(q);
	ZZ_pX phi;
	phi = to_ZZ_pX(modulus);
	ZZ_pE::init(phi);

	ZZ_pE f_, f_inv_;
	f_ = to_ZZ_pE(to_ZZ_pX(f));
	try{ f_inv_ = inv(f_); }
	catch(int e){ isfound = false; }
	ZZ_pX tp = rep(f_inv_);

	for(int i=0; i<degree; i++)
		SetCoeff(f_inv, i, rep(coeff(tp, i)));
}

void getClosestMod2(ZZ& out, ZZ& in, ZZ& p1, ZZ& p2)
{
	ZZ result;

	result = in  * p2 / p1;
	if(bit(result,0)!=bit(in,0))
	{

		if(result > 0)result ++;
		else result --;
	}
	out = result;
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



///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


