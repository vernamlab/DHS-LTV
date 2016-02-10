#ifndef PARAMS_H
#define PARAMS_H

#include "general.h"


class Params
{
    public:
                            Params();
            virtual        ~Params(){};

            void            AutoSetup();

            // GETs
            SetupType       setup_type()const{return setup_type_;};
            HomType         hom_type()const{return hom_type_;};
            RingType        ring_type()const{return ring_type_;};
            Flag            batch_flag()const{return batch_flag_;};
            Flag            fft_flag()const{return fft_flag_;};
            ReductionType   reduc_type()const{return reduc_type_;};
            int             pi()const{return pi_;};
            int             mu()const{return mu_;};
            int             nu()const{return nu_;};
            int             omega()const{return omega_;};
            int             sigma()const{return sigma_;};
            int             kappa()const{return kappa_;};
            int             lambda()const{return lambda_;};
            int             beta()const{return beta_;};
            int             depth()const{return depth_;};
            int             tau()const{return tau_;};
            int             rho()const{return rho_;};
            int             theta()const{return theta_;};
            double          delta()const{return delta_;};
            double          hermite()const{return hermite_;};



            // SETs
            void            set_setup_type(SetupType type);
            void            set_hom_type(HomType type);
            void            set_ring_type(RingType type);
            void            set_batch_flag(Flag flag);
            void            set_fft_flag(Flag flag);
            void            set_reduc_type(ReductionType type);
            void            set_pi(int val);
            void            set_mu(int val);
            void            set_nu(int val);
            void            set_omega(int val);
            void            set_sigma(int val);
            void            set_kappa(int val);
            void            set_lambda(int val);
            void            set_beta(int val);
            void            set_depth(int val);
            void            set_tau(int val);
            void            set_rho(int val);
            void            set_theta(int val);
            void            set_delta(double val);
            void            set_hermite(double val);


            // NOISE AND SECURITY ANALYSIS
            double          ComputeHermite();
            bool            CheckSecurity();
            bool            CheckNoiseGrowth();
            //void            ComputeInitialNoise();
            void            ComputeNoise(double &noise, double nu, double kappa, double tau, double omega, double adds, double err);
    protected:
    private:
            SetupType       setup_type_;
            HomType         hom_type_;
            RingType        ring_type_;
            Flag            batch_flag_, fft_flag_;
            ReductionType   reduc_type_;
            double          delta_, hermite_;
            int             pi_, mu_, nu_, sigma_, omega_, kappa_, beta_, lambda_, tau_, depth_, theta_, rho_;
};

#endif // PARAMS_H
