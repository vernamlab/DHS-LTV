#include "homlib.h"

Crypto::Crypto()
{
    gp_ = NULL;
    //fft_ek_ = NULL;
    srand(time(NULL));
	SetSeed(to_ZZ(time(NULL)));

homer_.Start();
    ParamSetup();
homer_.Stop();
homer_.ShowTime("Param Setup");

homer_.Start();
    RingSetup();
homer_.Stop();
homer_.ShowTime("Ring Setup");

homer_.Start();
    KeySetup();
homer_.Stop();
homer_.ShowTime("Key Setup");
}
Crypto::~Crypto()
{
    if(gp_ != NULL)
        delete gp_;
    /*if(fft_ek_.level_ != NULL)
    {
        for(int i=0; i<fft_ek_.length_; i++)
        {
            if(fft_ek_.level_[i].rep_list_ != NULL)
            {
                for(int j=0; j<fft_ek_.level_[i].block_cnt_; j++)
                    delete[] fft_ek_.level_[i].rep_list_[j];
                fft_ek_.level_[i].rep_list_;
            }
        }
        delete[] fft_ek_.level_;
    }*/

}

void Crypto::ParamSetup()
{
    // Setup gp_
    gp_ = new GlobalParams();
}
void Crypto::RingSetup()
{
    // Setup r_

    // 1. Poly Modulus, F(x)
    if(gp_->ring_type() == cyclotomic)
        CyclotomicPoly(r_.poly_, gp_->m());
    else
        XN_1(r_.poly_, gp_->n());

    // 2. Coeff modulus for each level, q_i
    r_.q_.SetLength(gp_->l()+1);

	ZZ temp;
	GenPrime(temp, gp_->i());
	while(temp%gp_->p() != 1)
        GenPrime(temp, gp_->i());

	r_.q_[gp_->l()] = temp;//Smallest q, last one

    if(gp_->hom_type() == she)
    {
        GenPrime(temp, gp_->k());
        while(temp%gp_->p() != 1)
            GenPrime(temp, gp_->k());
    }

	for(int i=gp_->l()-1; i>=0; i--)
		r_.q_[i] = (r_.q_[i+1])*temp;

    set_kappa(temp);

    // Helper polynomials for fast polynomial reduction
    // Set reduction_helpers_ for yarkin and barrett reduction
    // There is no helper poly for ntl reduction
    if(gp_->reduc_type() == yarkin)
    {
        int cnt = 0;
        ZZX mod = GetPolyModulus();
        int n = deg(mod);
        int temp = n;
        while(temp > 128)
        {
            cnt++;
            temp /= 2;
        }

        reduction_helpers_.SetLength(cnt);

        for(int i=0; i<cnt; i++)
        {
            SetCoeff(reduction_helpers_[i], n+(n>>(i+1)), 1);
            reduction_helpers_[i] = reduction_helpers_[i] % mod;
        }
    }
    else if(gp_->reduc_type() == barrett)
    {
        reduction_helpers_.SetLength(1);
        ZZX mod = GetPolyModulus();
        int n = deg(mod);
        SetCoeff(reduction_helpers_[0], 2*n-1, 1);
        reduction_helpers_[0] /= mod;
    }
}
void Crypto::KeySetup()
{
    if(gp_->hom_type() == fhe)
    {
        f_.SetLength(1);
        int cnt = (gp_->k()*(gp_->l()+1));
        cnt += gp_->w()-1;
        cnt /= gp_->w();
        ek_.SetLength(cnt);

        /*
        fft_ek_.length_ = gp_->l();
        fft_ek_.level_ = new FFTReps[fft_ek_.length_];
        */
    }
    else
    {
        ZZ total;
        power(total, 2, gp_->l());
        f_.SetLength(to_int(total));
    }

    ZZX f_inv;
    do
    {
        FindSecretKey();
    }while(!FindInverse(f_inv, f_[0], r_.poly_, r_.q_[0]));

    FindPublicKey(f_inv);

    if(gp_->hom_type() == she)
    {
        for(int i=1; i<f_.length(); i++)
        {
            f_[i] = f_[i-1] * f_[0];
            PolyCoeffReduce(f_[i], f_[i], gp_->ring_type(), r_.poly_, r_.q_[0]);
        }
    }
    else
    {
        ZZ temp = to_ZZ(1);
        temp = temp << gp_->w();    // 2^w
        ZZX msg = f_[0];
        for(int i=0; i<ek_.length(); i++)
        {
            encrypt(ek_[i], msg, 0);
            msg = msg * temp;
            CoeffReduce(msg, msg, r_.q_[0]);
        }
        //ComputeFFTek(0);
        //for(int i=1; i<fft_ek_.length_; i++)
        //{
            //Reduce EK
        //    ComputeFFTek(i);
        //}
    }
}

void Crypto::ComputeFFTek(int level)
{
    int cnt = ek_.length();             // leftblocks or word_blocksize or num_add
    int bits =  cnt * gp_->w();         // max_bitsize
    //cout << bits << endl;

    fft_ek_.level_[level].block_cnt_ = cnt;
    fft_ek_.level_[level].rep_list_ = new fftRep*[cnt];

    // How large of a prime do we need for an inner product?
    // Considering we will apply cnt number of mults
    // and cnt number of adds.
    long bound =  NumBits(gp_->n()+1) + bits + gp_->w() + NumBits(cnt) + 3;
    ZZ temp_prime;
    set(temp_prime);
    long prime_cnt = GetPrimeNumber(bound, temp_prime);
    zz_p::FFTInit(0);
    fft_ek_.level_[level].rep_cnt_ = prime_cnt;
    long k = NextPowerOfTwo(2*gp_->n()-1);

	for(int i=0; i<cnt; i++)
	{
        fft_ek_.level_[level].rep_list_[i] = new fftRep[prime_cnt];
        for(int j=0; j<fft_ek_.level_[level].rep_cnt_; j++)
        {
            zz_p::FFTInit(j);
            fft_ek_.level_[level].rep_list_[i][j].SetSize(k);
        }
        ComputeFFT(fft_ek_.level_[level].rep_list_[i], ek_[i], prime_cnt, gp_->n()-1);
    }
}

void Crypto::Encrypt(CipherText &out, const ZZX &in)
{
    Encrypt(out, in, 0);
}
void Crypto::Encrypt(CipherText &out, const ZZX &in, int index)
{
    out.set_fresh(true);
    out.set_level(index);
    out.set_multcnt(0);

    ZZX temp;
    encrypt(temp, in, index);
    out.set_ct(temp);
}
void Crypto::Decrypt(ZZX &out, CipherText &in)
{
    in.coeff_reduce(r_.q_[in.level()]);
    if(!in.fresh())
        Refresh(in);

    decrypt(out, in.ct(), in.d(), in.multcnt(), in.level());
}
void Crypto::ModSwitch(CipherText &ct)
{
homer_.Start();
    ZZX temp1, temp2;
    CoeffBalance(temp1, ct.ct(), r_.q_[ct.level()]);
    CoeffDivide(temp2, temp1, kappa_);
    CoeffParityMatch(temp2, temp1, gp_->p());

    ct.set_level(ct.level()+1);
    CoeffReduce(temp2, temp2, r_.q_[ct.level()]);
    ct.set_ct(temp2);

homer_.Stop();
homer_.ShowTime("ModSwitch");
}


void Crypto::Relin(CipherText &c)
{
homer_.Start();
    int cnt = ek_.length() - c.level()*(gp_->k()/gp_->w());
    vec_ZZX blocks;
    blocks.SetLength(cnt);
    SeparateBlocks(blocks, c.ct(), gp_->w());
    ZZX res;
    res = 0;
    for(int i=0; i<cnt; i++)
    {
        blocks[i] = blocks[i] * ek_[i];
        res += blocks[i];
    }
    PolyReduce(res, res, gp_->ring_type(), r_.poly_);
    CoeffReduce(res, res, r_.q_[c.level()]);
    c.set_ct(res);
homer_.Stop();
homer_.ShowTime("Relin");
}


/*
void Crypto::Relin(CipherText &c)
{
homer_.Start();

    zz_pBak bak;
	bak.save();
	int prime_cnt = fft_ek_.level_[c.level()].rep_cnt_;
    int cnt = fft_ek_.level_[c.level()].block_cnt_;
    long k = NextPowerOfTwo(2*gp_->n()-1);

    fftRep *temp_rep = new fftRep[prime_cnt];
    fftRep *out_rep = new fftRep[prime_cnt];
    FFTRepClear(temp_rep, prime_cnt);
    FFTRepClear(out_rep, prime_cnt);
    FFTRepSetK(temp_rep, prime_cnt, k);
    FFTRepSetK(out_rep, prime_cnt, k);

    vec_ZZX blocks;
    blocks.SetLength(cnt);
    SeparateBlocks(blocks, c.ct(), gp_->w());

    for(int i=0; i<cnt; i++)
    {
        ComputeFFT(temp_rep, blocks[i], prime_cnt, gp_->n()-1);
        for(int j=0; j<prime_cnt; j++)
        {
            zz_p::FFTInit(j);
            mul(temp_rep[j], temp_rep[j], fft_ek_.level_[c.level()].rep_list_[i][j]);
            add(out_rep[j], out_rep[j], temp_rep[j]);
        }
	}

    vec_zz_pX result_crt;
    result_crt.SetLength(prime_cnt);
    for(int i=0; i<prime_cnt; i++)
    {
		zz_p::FFTInit(i);
		FromfftRep(result_crt[i], out_rep[i], 0, 2*gp_->n()-2);
	}
	ZZX out;
	FFTtoPoly(out, result_crt, prime_cnt, gp_->n());
	CoeffReduce(out, out, r_.q_[c.level()]);

	delete[] temp_rep;
	delete[] out_rep;

	bak.restore();


homer_.Stop();
homer_.ShowTime("Relin");
}
*/

void Crypto::Refresh(CipherText &ct)
{
homer_.Start();
    //ct.poly_reduce(gp_->ring_type(), r_.poly_);
    PolyReduction(ct);
    if(gp_->hom_type() == fhe)
        Relin(ct);
    ModSwitch(ct);
    ct.set_fresh(true);
homer_.Stop();
homer_.ShowTime("Refresh");
}

void Crypto::PolyReduction(ZZX &out, const ZZX &in, ReductionType type)
{
    if(type == yarkin)
    {
        YarkinReduction(out, in, GetPolyModulus(), GetReductionHelpers());
    }
    else if(type == barrett)
    {
        ZZX u = reduction_helpers_[0];
        BarrettReduction(out, in, GetPolyModulus(), GetReductionHelpers()[0]);
    }
    else
    {
        out = in % GetPolyModulus();
    }

    out.normalize();
}

void Crypto::PolyReduction(CipherText &ct)
{
    ZZX temp = ct.ct();
    PolyReduction(temp, temp, gp_->reduc_type());
    ct.set_ct(temp);
}

void Crypto::Add(CipherText &out, CipherText &in1, CipherText &in2)
{
    if(!in1.fresh())
    {
        if(in2.fresh())
            //in1.poly_reduce(gp_->ring_type(), r_.poly_);
            PolyReduction(in1);
    }
    else if(!in2.fresh())
        //in2.poly_reduce(gp_->ring_type(), r_.poly_);
        PolyReduction(in2);

    out = in1 + in2;

    if(in2.level() > out.level())
        out.set_level(in2.level());
    out.coeff_reduce(r_.q_[out.level()]);

    out.set_fresh(in1.fresh()&&in2.fresh());

}


void Crypto::Mult(CipherText &out, CipherText &in1, CipherText &in2)
{
    if(!in1.fresh())
        Refresh(in1);
    if(!in2.fresh())
        Refresh(in2);

homer_.Start();
    out.set_level(in2.level());
    if(in1.level() > out.level())
        out.set_level(in1.level());


    ZZX temp = in1.ct() * in2.ct();
    PolyCoeffReduce(temp, temp, gp_->ring_type(), r_.poly_, r_.q_[out.level()]);
    CoeffReduce(temp, temp, r_.q_[out.level()]);
    out.set_ct(temp);

    if(gp_->hom_type() == fhe)
        out.set_multcnt(0);
    else
        out.set_multcnt(in1.multcnt() + in2.multcnt() + 1);

    out.set_fresh(false);
    out.set_d(in1.d()*in2.d());

homer_.Stop();
homer_.ShowTime("Mult Time");
}

void Crypto::Div(CipherText &out, const CipherText &in, const ZZ &d)
{
    int level = in.level();
    ZZ temp = InvMod(d, r_.q_[level]);
    out = in * temp;
    out.coeff_reduce(r_.q_[level]);
    out.set_d(out.d()*d);
}


CipherText& Crypto::Mult(CipherText &in1, CipherText &in2)
{
    Mult(in1, in1, in2);
    return in1;
}
CipherText& Crypto::Add(CipherText &in1, CipherText &in2)
{
    Add(in1, in1, in2);
    return in1;
}
CipherText& Crypto::Div(CipherText &in, const ZZ &d)
{
    Div(in, in, d);
    return in;
}


///////////////
// CIRCUITS
//////////////

void Crypto::Circuit_ZeroTest(CipherText &out, const CipherText &in)
{

}
void Crypto::Circuit_EqualityCheck(CipherText &out, const CipherText &in1, const CipherText &in2)
{
    out = in1 - in2;
    out.coeff_reduce(r_.q_[out.level()]);
    Circuit_ZeroTest(out, out);
}
//////////////
//////////////

void Crypto::FindSecretKey()
{
    RandPolyBalanced(f_[0], gp_->n(), gp_->b());
    f_[0] *= gp_->p();
    f_[0] += 1;
    CoeffReduce(f_[0], f_[0], r_.q_[0]);
}
void Crypto::FindPublicKey(const ZZX &f_inv)
{
    RandPolyBalanced(h_, gp_->n(), gp_->b());
    h_ *= f_inv;
    h_ *= gp_->p();
    PolyCoeffReduce(h_, h_, gp_->ring_type(), r_.poly_, r_.q_[0]);
}

void Crypto::encrypt(ZZX &out, const ZZX &in, int index)
{
homer_.Start();
    ZZX s, e;
    RandPolyBalanced(s, gp_->n(), gp_->b());
    RandPolyBalanced(e, gp_->n(), gp_->b());
    out = in;
    out = out + h_*s + e*gp_->p();
    PolyCoeffReduce(out, out, gp_->ring_type(), r_.poly_, r_.q_[index]);
homer_.Stop();
homer_.ShowTime("Encrypt");
}
void Crypto::decrypt(ZZX &out, const ZZX &in, const ZZ &divisor, int index, int level)
{
homer_.Start();
    out = in * f_[index];
    PolyReduce(out, out, gp_->ring_type(), r_.poly_);
    CoeffBalance(out, out, r_.q_[level]);
    CoeffReduce(out, out, to_ZZ(gp_->p())/divisor);
homer_.Stop();
homer_.ShowTime("Decrypt");
}
