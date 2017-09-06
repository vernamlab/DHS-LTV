#include "ltv.h"

//////////////////////////////////////////////////////////////////////////////////////////////
ltv::ltv(int dm, int num, int max_bit, int tableS, ZZ qq, int wordSize){

	word 	= wordSize;
	ZZ BB 	= to_ZZ("1");
	Set(gp, qq, BB, dm);

	p 			= gp.q;
	B 			= gp.B;
	N 			= gp.N_PolyDegree;
	degree_m 	= gp.M_CycDegree;
	rand_seed 	 = 0;
	message_rand = 0;
	reducsetflag = false;

	clear(modulus, N+1);
	ZZ_p::init(p);
	logq = max_bit/num;
/*
	init_bitsize = max_bitsize;
	tableSize 	 = tableS;

	max_bitsize = max_bit;

	pk = new ZZX[num];
	sk  = new ZZX[num];
	q_list = new ZZ[num];
	ek2 = new struct eval_key[num];

	for(int i=0; i<num; i++)
		ek2[i].key 	= new ZZX[max_bitsize];


	fftKeys = new FFTPolyList[tableSize];
	for(int i=0; i<tableSize; i++){
		fftKeys[i].polynum		= max_bitsize;
		fftKeys[i].polyList 	= new PolyFFTRep[max_bitsize];
		for(int j=0; j<max_bitsize; j++)
			fftKeys[i].polyList[j].poly_rep = new fftRep[5];
	}
*/
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
	logq = max_bit/num;
/*
	init_bitsize = max_bit;
	tableSize = tableS;

	max_bitsize = max_bit;

	pk = new ZZX[num];
	sk  = new ZZX[num];
	q_list = new ZZ[num];
	ek2 = new struct eval_key[num];

	for(int i=0; i<num; i++)
		ek2[i].key 	= new ZZX[max_bitsize];

	fftKeys = new FFTPolyList[tableSize];
	for(int i=0; i<tableSize; i++){
		fftKeys[i].polynum		= max_bitsize;
		fftKeys[i].polyList 	= new PolyFFTRep[max_bitsize];
		for(int j=0; j<max_bitsize; j++)
			fftKeys[i].polyList[j].poly_rep = new fftRep[5];
	}
*/

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

	//Arith_PolyReduce(result, result);
	result = Arith_PolyCoeffReduce(result, i);
	//coeff_reduction(result, result, q_list[i], N);

	return result;
}

ZZX ltv::Prim_Decrypt(ZZX c, int index){

	ZZ x;
	ZZX result, t;

	t = sk[index]*c;
	t = Arith_PolyCoeffReduce(t, index);
//	coeff_reduction(t, t, q_list[index], N);
	for(int i=0; i<N; i++){
		x = coeff(t, i);
		if(x>((q_list[index]-1)/2))
			x = x-q_list[index];
		SetCoeff(result, i, (x%p));
	}
	return result;
}

ZZX ltv::Prim_Decrypt2(ZZX c, int index){

	ZZ x;
	ZZX result, t;

	t = sk[index]*sk[index]*c;
	t = t%modulus;
	t = Arith_PolyCoeffReduce(t, index);
//	coeff_reduction(t, t, q_list[index], N);
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
	int leftWords = max_bitsize-index*logq;
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

	ZZX modin, modout;
	for(int i=0; i<N; i++){
		tp = coeff(x, i);
		if(tp>((q_list[index]-1)/2))
			tp = tp - q_list[index];

		SetCoeff(modin, i, tp%p);
		RoundModP(tp, tp, q_list[index], q_list[index+1], p);
		SetCoeff(modout, i, tp%p);

		SetCoeff(res, i, tp%q_list[index+1]);
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

	int leftBits = max_bitsize-level*logq;
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
		if(rand()%2)SetCoeff(result,i, tp);
		else SetCoeff(result,i, 0-tp);
	}
	return result;
}

void ltv::Func_ModulusFindRing(int num, int max_bit, int diff, ZZX modu){

	ZZ 	t, co, k;
	max_bitsize = max_bit;
	init_bitsize = max_bitsize;

	GenPrime(t, diff);
	k = t/p;
	t = k*p+1;
	while(ProbPrime(t, 80) != 1){
		k++;
		t = k*p+1;
	}
	co = t;
	for(int i=0; i<num; i++){
		q_list[num-1-i] = t;
		t=t*co;
	}
//	for(int i=0; i<num; i++)
//			cout << q_list[i] << endl;

	/*
	ZZ q[num];
	int bs;
	bs = diff;
	for(int i=0; i<num; i++) {
		ZZ 	t, c;
		max_bitsize = max_bit;
		init_bitsize = max_bitsize;
		GenPrime(t, bs);
		c = t/p;
		t = c*p+1;
		while(ProbPrime(t) != 1){
			c++;
			t = c*p+1;
		}
		q[i] = t;
		bs = bs + diff;
	}

	for(int i=0; i<num; i++)
		q_list[i] = q[num-1-i];
	
	for(int i=0; i<num; i++) 
		cout << q_list[i] << endl;
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
#if PolyModXN_1 == 1
	SetCoeff(modulus, N, 1);
	SetCoeff(modulus, 0, -1);
#else
	modulus = ComputeFastCycModulus(degree_m);
#endif
}

void ltv::Func_ComputeKeysRingRelinFFT(int num, int tblsize){
	tableSize = tblsize;

	if(reducsetflag == false)
		SetFastPolyModulus();

////////////////////////First Key//////////////////////////////////////
	ZZX f, g, ft, f_inv;
	bool isfound = false;


	ft  = Func_Sample();
	ZZX ft_ = ft;
	int bitsize = max_bitsize;
//	for(int i=0; i<num; i++){
		int i=0;
		while(!isfound){
			isfound = true;
			coeff_reduction(ft, ft_, q_list[i], N);
			f = ft*p + 1;
			coeff_reduction(f, f, q_list[i], N);
			find_inverse(f_inv, f, q_list[i], N, isfound, modulus);
			coeff_reduction(f_inv, f_inv, q_list[i], N);
		}
		key_f = f;
		key_fi = f_inv;

		ZZX eee = f_inv * f;
		Arith_PolyReduce(eee, eee);
		coeff_reduction(eee, eee, q_list[i], N);

		isfound = false;
		g 	  = Func_Sample();
		coeff_reduction(g, g, q_list[i], N);
		pk[i] = g*f_inv;
		Arith_PolyReduce(pk[i], pk[i]);

		pk[i] = pk[i]*p;
		sk[i] = f;


		coeff_reduction(pk[i], pk[i], q_list[i], N);
		coeff_reduction(sk[i], sk[i], q_list[i], N);

		ZZX	tppk = pk[i];
		ZZX tpek = sk[i];
		compute_eval_onerelin(ek2[i], tppk, tpek, i);
	/////////////////////////
		ConvertKeystoFFT(ek2[i], fftKeys[i], max_bitsize, word_blocksize);
		bitsize = bitsize - logq;


		for(i=1; i<num; i++){
			coeff_reduction(pk[i], pk[i-1], q_list[i], N);
			coeff_reduction(sk[i], sk[i-1], q_list[i], N);
		}

		for(i=0; i<num; i++){
		ZZX keys_ = pk[i]*sk[i];
		ZZX t;
		t = Arith_PolyCoeffReduce(keys_, i);

		for(int j=0; j<N; j++){
			ZZ tp = coeff(t, j);
			if(tp>((q_list[i]-1)/2))
				tp = tp - q_list[i];
				SetCoeff(t, j, tp);
			}
			coeff_reduction(t, t, p, N);
		}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void PolyRedXN_1(ZZX &out, ZZX &in, int N){
	ZZX res;
	for(int i=0; i<N; i++){
		SetCoeff(res, i, coeff(in,i)+coeff(in,i+N));
	}
	out = res;
}

void ltv::Arith_PolyReduce(ZZX &out, ZZX &in){
#if PolyModXN_1 == 0
	myr.Reduce(out, in);
#elif PolyModXN_1 == 1
	PolyRedXN_1(out, in, N);
#elif PolyModXN_1 == 2
	PolyRedXN_1(out, in, N+1);
	for(int i=0; i<N+1; i++)
		SetCoeff(out, i, coeff(out,i)-coeff(out, N));
#endif
}

void ltv::Arith_CoeffReduce(ZZX &out, ZZX &in, int index){
	int l = in.rep.length();
	coeff_reduction(out, in, q_list[index], l);
}

ZZX ltv::Arith_PolyCoeffReduce(ZZX &in, int index){
	ZZX r;
	r = in;
	Arith_CoeffReduce(r, r, index);
	Arith_PolyReduce(r, r);
	Arith_CoeffReduce(r, r, index);
	return r;
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


void ltv::IO_ReadModulus(fstream &fs){
	fs >> modulus;
	if(reducsetflag == false)
		SetFastPolyModulus();
}

void ltv::IO_ReadQs(fstream &fs, int s){
	for(int i=0; i<s; i++)
		fs >> q_list[i];
}

void ltv::IO_ReadKeys(fstream &fs, int s){
	fs >> pk[0];
	for(int i=0; i<s; i++)
		fs >> sk[i];

	for(int i=0; i<s; i++)
		fs >> ek2[0].key[i];

	if(reducsetflag == false)
		SetFastPolyModulus();

	ConvertKeystoFFT(ek2[0], fftKeys[0], max_bitsize, max_bitsize);
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

void coeff_reduction_q_2(ZZX &out, ZZX &in, ZZ &q, int N){
	ZZ t;
	ZZ s = (q-1)/2;
	for(int i=0; i<N; i++){
		t = coeff(in, i);
		if(t>s)
			t = t-q;
		SetCoeff(out, i, t);
	}
}

void getClosestMod2(ZZ& out, ZZ& in, ZZ& p1, ZZ& p2)
{
	ZZ result;

	result = in  * p2 / p1;
	if(bit(result,0)!=bit(in,0)){

		if(result > 0)result ++;
		else result --;
	}
	out = result;
}

void RoundModP(ZZ& out, ZZ& in, ZZ& p1, ZZ& p2, ZZ& modp)
{
	ZZ result, tm0, tm1;

	result = in  * p2 / p1;


	tm0 = result%modp;
	tm1 = in%modp;

	if(result > 0){
		if(tm0>tm1)
			result = result + tm1 - tm0;
		else if(tm0<tm1)
			result = result - tm0 + tm1;
	}
	else{
		if(tm0>tm1)
			result = result + tm1 - tm0;
		else if(tm0<tm1)
			result = result - tm0 + tm1;
	}
//	cout << result << endl;
	out = result;
}

ZZ ComputePR(ZZ q, ZZ m){
	ZZ res, l;
	ZZ t,c;

	t = (q-1)/m;
	l = to_ZZ("2");

	bool is  = true;

	while(is){
		c = PowerMod(l, t, q);

		if(c == to_ZZ("1"))
			l++;
		else{
			is = false;
			vec_ZZ f = findFactors(m);
			for(int i=0; i<f.length(); i++){
				if(PowerMod(c, m/f[i], q) == to_ZZ("1")){
					l++;
					is = true;
				}
			}
		}
	}
	res = l;
	res = PowerMod(res, (q-1)/m, q);
	return res;
}

vec_ZZ findFactors(ZZ f){

	ZZX v;
	vec_ZZ res;

	ZZ s = to_ZZ("2");

	int l = 0;
	while(s<f){
		if(GCD(s,f) != to_ZZ("1")){
			SetCoeff(v, l, s);
			l++;
		}
		s++;
		s = NextPrime(s);
	}

	res.SetLength(l);
	for(int j=0; j<l; j++)
		res[j] = coeff(v, j);

	return res;
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




/*
void ltv::ModulusFind(int num, int max_bit, int diff){

	ZZ 	t, co;
	int bit_size = max_bit;
	max_bitsize = max_bit;
	init_bitsize = max_bitsize;

	q_list = new ZZ[num];
	for(int i=0; i<num; i++){
		t  = NextPrime(power2_ZZ(bit_size));
		co = t/(p*degree_m);//co = t/(p*N);
		bool loop = true;
		while(loop){
			t = co*(p*degree_m)+1;//t = co*(p*N)+1;
			if(ProbPrime(t, 100)){
				loop = false;
				if(GCD(t-1, to_ZZ(degree_m)) == to_ZZ(degree_m))
					loop = false;
				else
					loop = true;
			}
			else
				co++;
		}
		bit_size = bit_size - diff;
		q_list[i] = t;
	}

//	pr_list	= new ZZ[num];

//	for(int i=0; i<num; i++)
//		pr_list[i] = ComputePR(q_list[i], to_ZZ(degree_m));
}
*/
/*
void ltv::ComputeKeys(int num){

	pk = new ZZX[num];		sk  = new ZZX[num];
	ek2 = new struct eval_key[num-1];

	for(int i=0; i<num-1; i++){
		ek2[i].key 	= new ZZX[max_bitsize];
	}

	if(reducsetflag == false)
		SetFastPolyModulus();

    ///////////////////First Key///////////////////////
	ZZX f, g, ft, f_inv;
	bool isfound = false;

	while(!isfound){
		isfound = true;
		ft  = Sample();
		coeff_reduction(ft, ft, q_list[0], N);

		f = ft*p + 1;
		coeff_reduction(f, f, q_list[0], N);
		find_inverse(f_inv, f, q_list[0], N, isfound, modulus);
		coeff_reduction(f_inv, f_inv, q_list[0], N);
	}
	isfound = false;

	g 	  = Sample();
	coeff_reduction(g, g, q_list[0], N);
	pk[0] = MulMod(g, f_inv, modulus);

	pk[0] = pk[0]*p;
	sk[0] = f;

	coeff_reduction(pk[0], pk[0], q_list[0], N);
	coeff_reduction(sk[0], sk[0], q_list[0], N);

#ifdef PrintKeysOut
	cout << "Keys 0 Complete!" << endl;
#endif

    //////////////////////Rest Key//////////////////////
	for(int i=1; i<num; i++){
		while(!isfound){
			isfound = true;
			ft  = Sample();
			coeff_reduction(ft, ft, q_list[i-1], N);

			f = ft*p + 1;
			coeff_reduction(f, f, q_list[i-1], N);

			find_inverse(f_inv, f, q_list[i-1], N, isfound, modulus);
			coeff_reduction(f_inv, f_inv, q_list[i-1], N);
		}
		isfound = false;

		g 	  = Sample();
		coeff_reduction(g, g, q_list[i-1], N);

		ZZX tppk 	= MulMod(g, f_inv, modulus);
		tppk 		= tppk*p;
		sk[i]		= f;
		pk[i]		= tppk;

		coeff_reduction(tppk, tppk, q_list[i-1], N);
		coeff_reduction(sk[i], sk[i], q_list[i-1], N);

		coeff_reduction_q_2(sk[i], sk[i], q_list[i-1], N);

		ZZX tpek = sk[i-1];
		compute_eval(ek2[i-1], tppk, tpek, i-1);

		#ifdef PrintKeysOut
			cout << "Keys " << i << " Complete!" << endl;
		#endif
	}
}
*/

/*
void ltv::ComputeKeysRingNoRelin(int num){
	pk = new ZZX[num];		sk  = new ZZX[num];

	if(reducsetflag == false)
		SetFastPolyModulus();

    ////////////////////First Key///////////////////////////
	ZZX f, g, ft, f_inv;
	bool isfound = false;

	while(!isfound){
		isfound = true;
		ft  = Sample();
		coeff_reduction(ft, ft, q_list[0], N);

		f = ft*p + 1;
		coeff_reduction(f, f, q_list[0], N);
		find_inverse(f_inv, f, q_list[0], N, isfound, modulus);
		coeff_reduction(f_inv, f_inv, q_list[0], N);
	}
	isfound = false;

	g 	  = Sample();
	coeff_reduction(g, g, q_list[0], N);
	pk[0] = MulMod(g, f_inv, modulus);

	pk[0] = pk[0]*p;
	sk[0] = f;

	coeff_reduction(pk[0], pk[0], q_list[0], N);
	coeff_reduction(sk[0], sk[0], q_list[0], N);

#ifdef PrintKeys
	cout << "Keys 0 Complete!" << endl;
#endif

    /////////////////////Rest Key/////////////////////////
	for(int i=1; i<num; i++){
		sk[i] = (sk[i-1]*sk[i-1]);
		PolyReduce(sk[i], sk[i]);
		coeff_reduction(sk[i], sk[i], q_list[i], N);

		#ifdef PrintKeys
			cout << "Keys " << i << " Complete!" << endl;
		#endif
	}
}
*/
/*
ZZX ltv::Relin(ZZX c, int index){
	int bits[to_long(N)][max_bitsize];
	ZZ t;
	for(int i=0; i<N; i++){
		t = coeff(c,i);
		for(int j=0; j<max_bitsize; j++){
			bits[i][j] = bit(t, j);
		}
	}

	ZZ tp[2*to_long(N)];

	for(int i=0; i<2*to_long(N); i++)
		tp[i] = to_ZZ("0");

	for(int b=0; b<max_bitsize; b++){
		for(int i=0; i<N; i++){
			if(bits[i][b] == 1){
				for(int j=0; j<N; j++)
					tp[j+i] = tp[j+i] + coeff(ek2[index].key[b], j);
			}
		}
	}

	ZZX result;
	for(int i=0; i<2*N; i++)
		SetCoeff(result, i, tp[i]);

	result = result%modulus;
	for(int i=0; i<N; i++)
		coeff_reduction(result, result, q_list[index], N);

	return result;
}
*/

/*
void ltv::ComputeKeysRingRelin(int num){
	pk = new ZZX[num];		sk  = new ZZX[num];
	ek2 = new struct eval_key[num];

	for(int i=0; i<num; i++){
		ek2[i].key 	= new ZZX[max_bitsize];
	}

	if(reducsetflag == false)
		SetFastPolyModulus();

////////////////////////First Key//////////////////////////////////////
	ZZX f, g, ft, f_inv;
	bool isfound = false;

	while(!isfound){
		isfound = true;
		ft  = Sample();
		coeff_reduction(ft, ft, q_list[0], N);
		f = ft*p + 1;
		coeff_reduction(f, f, q_list[0], N);
		find_inverse(f_inv, f, q_list[0], N, isfound, modulus);
		coeff_reduction(f_inv, f_inv, q_list[0], N);
	}
	isfound = false;
	g 	  = Sample();
	coeff_reduction(g, g, q_list[0], N);
	pk[0] = MulMod(g, f_inv, modulus);

	pk[0] = pk[0]*p;
	sk[0] = f;

	coeff_reduction(pk[0], pk[0], q_list[0], N);
	coeff_reduction(sk[0], sk[0], q_list[0], N);

	ZZX	tppk = pk[0];
	ZZX tpek = sk[0];

	compute_eval_onerelin(ek2[0], tppk, tpek, 0);

	for(int j=1; j<num; j++){
		for(int i=0; i<max_bitsize; i++){
			coeff_reduction(ek2[j].key[i], ek2[0].key[i], q_list[j], N);
		}
	}

#ifdef PrintKeys
	cout << "Keys 0 Complete!" << endl;
#endif

//////////////////////////Rest Key///////////////////////////////////
	for(int i=1; i<num; i++){
		sk[i] = sk[i-1];
		coeff_reduction(sk[i], sk[i], q_list[i], N);

		#ifdef PrintKeys
			cout << "Keys " << i << " Complete!" << endl;
		#endif
	}
}
*/

/*
ZZX ntru::RelinRing(ZZX c, int index){
	ZZX bits;

	ZZX tp;
	for(int j=0; j<max_bitsize; j++){
		for(int i=0; i<N; i++)
			SetCoeff(bits, i, bit(coeff(c,i),j));

		tp = tp + bits*ek2[index].key[j];
	}

	int y = 2*N;
	coeff_reduction(tp, tp, q_list[index], y);

	return tp;
}
*/
/*
	pk = new ZZX[num];		sk  = new ZZX[num];		ek2 = new struct eval_key;

	ek2[0].key 	= new ZZX[max_bitsize];

	fftKeys = new FFTPolyList[tableSize];
	for(int i=0; i<tableSize; i++){
		fftKeys[i].polynum		= max_bitsize;
		fftKeys[i].polyList 	= new PolyFFTRep[max_bitsize];
	}
*/
/*
void ntru::RingBasedKeygen(int num){

	if(reducsetflag == false)
		SetFastPolyModulus();

	ZZX f, g, ft, f_inv;
	bool isfound = false;

	while(!isfound){
		isfound = true;
		ft  = Sample();
		coeff_reduction(ft, ft, q_list[0], N);

		f = ft*p + 1;
		coeff_reduction(f, f, q_list[0], N);
		find_inverse(f_inv, f, q_list[0], N, isfound, modulus);
		coeff_reduction(f_inv, f_inv, q_list[0], N);
	}
	isfound = false;

	g 	  = Sample();
	coeff_reduction(g, g, q_list[0], N);
	pk[0] = g*f_inv;
	PolyReduce(pk[0], pk[0]);
	pk[0] = pk[0]*p;
	sk[0] = f;

	coeff_reduction(pk[0], pk[0], q_list[0], N);
	coeff_reduction(sk[0], sk[0], q_list[0], N);

	for(int i=1; i<num; i++){
		sk[i] = (sk[i-1]*sk[i-1]);
		PolyReduce(sk[i], sk[i]);
		coeff_reduction(sk[i], sk[i], q_list[i], N);
	}
}
*/
/*
void ntru::RingBasedKeygen(int num){

	if(reducsetflag == false)
		SetFastPolyModulus();

	ZZX f, g, ft, f_inv;
	bool isfound = false;

	while(!isfound){
		isfound = true;
		ft  = Sample();
		coeff_reduction(ft, ft, q_list[0], N);

		f = ft*p + 1;
		coeff_reduction(f, f, q_list[0], N);
		find_inverse(f_inv, f, q_list[0], N, isfound, modulus);
		coeff_reduction(f_inv, f_inv, q_list[0], N);
	}
	isfound = false;

	g 	  = Sample();
	coeff_reduction(g, g, q_list[0], N);
	pk[0] = g*f_inv;
	PolyReduce(pk[0], pk[0]);
	pk[0] = pk[0]*p;
	sk[0] = f;

	coeff_reduction(pk[0], pk[0], q_list[0], N);
	coeff_reduction(sk[0], sk[0], q_list[0], N);

	for(int i=1; i<num; i++){
		sk[i] = (sk[i-1]*sk[i-1]);
		PolyReduce(sk[i], sk[i]);
		coeff_reduction(sk[i], sk[i], q_list[i], N);
	}
}
*/
/*
ZZX	ntru::FFTTestFunc(fftRep *R, int priNum){
	ZZX ffttest;
	zz_pX *rescrt;
	rescrt = new zz_pX[priNum];

	for(int i=0; i<priNum; i++){
		zz_p::FFTInit(i);
		FromfftRep(rescrt[i], R[i], 0, N+N-2);
	}
	int t = priNum;
	FromCrt2Poly(ffttest, rescrt, t);

	 return ffttest;
}
*/
/*
ZZX ntru::ModSwitchX_M_1(ZZX &x, int index){
	ZZX res;
	ZZ tp;
	for(int i=0; i<degree_m; i++){
		tp = coeff(x, i);
		if(tp>((q_list[index]-1)/2))
			tp = tp - q_list[index];

		getClosestMod2(tp, tp, q_list[index], q_list[index+1]);
		SetCoeff(res, i, tp);
	}
	return res;
}

ZZX	ntru::PolyModX_N_1(ZZX &in){
	ZZX res;
	for(int i=0; i<degree_m; i++){
		SetCoeff(res, i, coeff(in, i) - coeff(in, i+degree_m));
	}
	return res;
}

void ntru::compute_eval_XN1(eval_key &ek2, ZZX tppk, ZZX tpek, int index){

	ZZX tpek2 = tpek*tpek;
	tpek2 = PolyModX_N_1(tpek2);
	coeff_reduction(tpek2, tpek2, q_list[index], degree_m);

	ZZX s, e, result, tp;
	for(int i=0; i<max_bitsize; i++){

		tp = tpek2 * power2_ZZ(i);
		s = Sample();
		e = Sample();
		coeff_reduction(s, s, q_list[index], N);
		coeff_reduction(e, e, q_list[index], N);
		coeff_reduction(tp, tp, q_list[index], degree_m);

		ek2.key[i] = tppk*s + e*p + tp;
		ek2.key[i] = PolyModX_N_1(ek2.key[i]);

		coeff_reduction(ek2.key[i], ek2.key[i], q_list[index], degree_m);
	}
}

ZZX ntru::Encrypt_XN1(ZZX m, int i){

	ZZX s, e, result;
	s = Sample();
	e = Sample();
	coeff_reduction(s, s, q_list[i], N);
	coeff_reduction(e, e, q_list[i], N);

	result = pk[i]*s+ e*p + m;//MulMod(pk[i], s, modulus) + e*p + m;
	result = PolyModX_N_1(result);
	coeff_reduction(result, result, q_list[i], degree_m);

	return result;
}

ZZX ntru::MulModZZX_XN1(ZZX &a, ZZX &b, int index){
	ZZX t = a*b;//MulMod(a, b, modulus);
	t = PolyModX_N_1(t);
	coeff_reduction(t, t, q_list[index], degree_m);
	return t;
}

ZZX ntru::Relin_XN1(ZZX c, int index){
	int bits[to_long(degree_m)][max_bitsize];
	ZZ t;
	for(int i=0; i<degree_m; i++){
		t = coeff(c,i);
		for(int j=0; j<max_bitsize; j++){
			bits[i][j] = bit(t, j);
		}
	}

	ZZ tp[2*to_long(degree_m)];

	for(int i=0; i<2*to_long(degree_m); i++)
		tp[i] = to_ZZ("0");

	for(int b=0; b<max_bitsize; b++){
		for(int i=0; i<degree_m; i++){
			if(bits[i][b] == 1){
				for(int j=0; j<degree_m; j++)
					tp[j+i] = tp[j+i] + coeff(ek2[index].key[b], j);
			}
		}
	}

	ZZX result;
	for(int i=0; i<2*degree_m; i++)
		SetCoeff(result, i, tp[i]);

	result = PolyModX_N_1(result);
	coeff_reduction(result, result, q_list[index], degree_m);

	return result;
}


ZZX ntru::ModSwitch_XN1(ZZX &x, int index){

	ZZ tp;
	for(int i=0; i<degree_m; i++){
		tp = coeff(x, i);
		if(tp>((q_list[index]-1)/2))
			tp = tp - q_list[index];

		getClosestMod2(tp, tp, q_list[index], q_list[index+1]);
		SetCoeff(x, i, tp);
	}
	return x;
}


ZZX ntru::Decrypt_XN1(ZZX c, int index){

	ZZ x;
	ZZX result;
	ZZX t = sk[index]*c;//MulMod(sk[index], c, modulus);
	t = PolyModX_N_1(t);
	coeff_reduction(t, t, q_list[index], degree_m);

	for(int i=0; i<degree_m; i++){
		x = coeff(t, i);
		if(x>((q_list[index]-1)/2))
			x = x-q_list[index];
		SetCoeff(result, i, (x%p));
	}
	return result;
}

ZZX	ntru::PolyModX_Np1(ZZX &in){
	ZZX res;
	for(int i=0; i<degree_m; i++){
		SetCoeff(res, i, coeff(in, i) + coeff(in, i+degree_m));
	}
	return res;
}

ZZX	ntru::MulModPolyX_M_1(ZZX &a, ZZX &b, int index){
	ZZX t;
	mul(t, a, b);
	t = PolyModX_Np1(t);
	coeff_reduction(t, t, q_list[index], degree_m);
	return t;
}
*/
