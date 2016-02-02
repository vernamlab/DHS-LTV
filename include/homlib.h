#ifndef HOMLIB_H_
#define HOMLIB_H_

#include "general.h"

using namespace std;
using namespace NTL;


class Crypto{
public:

///////////////////////////////////////////////////////////////////////////////
        Crypto();
        ~Crypto();

        //void    ParamSetup();
        //void    KeySetup();

        const ZZX&      GetSecretKey(int index)const {return f_[index];};
        const ZZX&      GetPublicKey()const {return h_;};
        const ZZX&      GetPolyModulus()const {return r_.poly_;};
        const vec_ZZX&  GetReductionHelpers()const {return reduction_helpers_;};
        const ZZ&       GetCoeffModulus(int index)const {return r_.q_[index];};
        const ZZ&       GetMessageModulus()const {return to_ZZ(gp_->p());};
        RingType        GetRingType()const {return gp_->ring_type();};


        void    Encrypt(CipherText &out, const ZZX &in);
		void    Encrypt(CipherText &out, const ZZX &in, int index);
		void	Decrypt(ZZX &out, CipherText &in);
		void 	ModSwitch(CipherText &ct);
		void    Relin(CipherText &c);
		void    Refresh(CipherText &ct);

        void    Mult(CipherText &out, CipherText &in1, CipherText &in2);
        void    Add(CipherText &out, CipherText &in1, CipherText &in2);
        void    Div(CipherText &out, const CipherText &in, const ZZ &d);

        CipherText &    Mult(CipherText &in1, CipherText &in2);
        CipherText &    Add(CipherText &in1, CipherText &in2);
        CipherText &    Div(CipherText &in, const ZZ &d);

        const ZZ&   kappa()const{return kappa_;};
        void        set_kappa(const ZZ &k){kappa_ = k;};

        void    Circuit_ZeroTest(CipherText &out, const CipherText &in);
        void    Circuit_EqualityCheck(CipherText &out, const CipherText &in1, const CipherText &in2);

        void    TimerStart(){homer_.Start();};
        void    TimerStop(){homer_.Stop(); homer_.ShowTime("Time ");};

        void    PolyReduction(CipherText &out, const CipherText &in);
        void    PolyReduction(ZZX &out, const ZZX &in);




private:
        ZZX             h_;
        vec_ZZX         f_, ek_;		// list of public (h), secret (f), and evaluation (ek) keys for each level
        FFTEvalKey      fft_ek_;
        //fftRep          ek_fft_;
        GlobalParams     *gp_;
        Ring            r_;
        vec_ZZX         reduction_helpers_;
        myTimer         homer_;
        ZZ              kappa_;         // cutting factor


        void    ParamSetup();
        void    KeySetup();
        void    RingSetup();

        void    FindSecretKey();
        void    FindPublicKey(const ZZX &f_inv);

        void    ComputeFFTek(int level);

        void    encrypt(ZZX &out, const ZZX &in, int index);
		void	decrypt(ZZX &out, const ZZX &in, const ZZ &divisor, int index, int level);

		void    PolyReduction(ZZX &out, const ZZX &in, ReductionType type);

};

#endif /* HOMLIB_H_ */
