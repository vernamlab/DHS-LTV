#include "params.h"


Params::Params()
{
    // APPLICATION
    #ifdef PI
        int p = pow(2, PI);
        set_pi(p);                                     // Message Space Z_{p}
    #endif // PI

    #ifdef CIRCUIT_DEPTH
        set_depth(CIRCUIT_DEPTH);                       // Obvious?
    #endif // CIRCUIT_DEPTH

    #ifdef BETA
        int b = pow(2, BETA);
        set_beta(b);                                  // Initial noise bound B
    #endif // BETA

    // HOMO TYPE
    #ifdef SHE
        set_hom_type(she);
    #elif defined LHE
        set_hom_type(lhe);
    #elif defined FHE
        set_hom_type(fhe);
    #endif // FHE

    // BATCHED?
    #ifdef BATCH
        set_batch_flag(on);
    #else
        set_batch_flag(off);
    #endif // BATCH

    // FFT speed for relinearization?
    #ifdef FFT
        set_fft_flag(on);
    #else
        set_fft_flag(off);
    #endif // FFT

    // Polynomial Reduction Algorithm
    #ifdef FAST_REDUCTION
        set_reduc_type(fast);
    #elif defined BARRETT_REDUCTION
        set_reduc_type(barrett);
    #elif defined NTL_REDUCTION
        set_reduc_type(ntl);
    #endif // NTL_REDUCTION

    // Ring Type
    #ifdef CYCLOTOMIC
        set_ring_type(cyclotomic);
    #elif defined XN_1
        set_ring_type(xn_1);
    #endif // CYCLOTOMIC

    // SETUP
    #ifdef MANUAL_SETUP
        set_setup_type(manual);

        // NOISE
        set_kappa(KAPPA);                               // Cutting factor bit size K
        set_lambda(LAMBDA);                             // Final step (min) coefficient q size
        set_sigma(SIGMA);                               // First step (max) coefficient q size

        #ifdef OMEGA
            set_omega(OMEGA);                           // Relin block bit size
            set_tau(TAU);                               // Eval Keys will be : E(f*2^{OMEGA_*i}) where i=0,...,TAU_-1
        #else
            set_omega(kappa()/(2*(depth()+1)));         // If not defined in def.h, this is the optimized relin block size
            set_tau((sigma()+omega()-1)/omega());
        #endif // OMEGA

        #ifdef CYCLOTOMIC
            #ifndef BATCH
                set_mu(MU);
                set_nu(EulerToient(MU));
            #else
                #if PI == 1
                    set_mu(MU);
                    set_nu(EulerToient(MU));
                    set_theta(ComputeFactorDegree(MU)); // T : Factor Degree. i.e. F(x) can be factored into a product of equal degree polynomials
                #else
                    set_theta(2);
                    set_mu(p()*p()-1);
                    set_nu(EulerToient(MU));
                #endif // PI
            #endif // BATCH
        #else   //XN_1
            set_nu(NU);
            #ifdef BATCH
                set_theta(1);                           // T : Factor Degree. i.e. F(x) can be factored into a product of equal degree (13) polynomials
                set_pi(FindPrimeCongOne(NU));           // p = 1 mod n
            #endif // BATCH
        #endif // CYCLOTOMIC

        #ifdef BATCH
            set_rho(nu()/theta());                    // R : Factor (Message Slot) Count.
        #else
            set_rho(1);
        #endif // BATCH

        #ifndef SIGMA
            set_sigma(lambda() + kappa() * depth());    // Largest (first level) q bit size for the ring R = Z_{q_{0}}/F(x)
                                                        // SIGMA_ = log_2(q_{0})
        #endif // SIGMA

        #ifdef SECURE
            set_delta(DELTA);
            CheckSecurity();
        #endif // SECURE
    #elif defined AUTO_SETUP
        set_setup_type(automatic);
        AutoSetup();
    #endif // MANUAL_SETUP

}

void Params::AutoSetup()
{
    #ifdef DEBUG_PARAM_SETUP
        cout << endl << endl << "AUTO PARAM SETUP START" << endl<< endl;
    #endif // DEBUG_PARAM_SETUP

    #ifdef SECURE
        set_delta(DELTA);
        int kappa = 20;
        int sigma = kappa * (depth()+1);
        int nu = pow(2,10);

        bool found = false;

        while(!found)
        {
            while(sigma > (log2(delta_)*4*nu) )
                nu *= 2;
            #ifdef DEBUG_PARAM_SETUP
                cout << "Secure n : " << nu << endl;
            #endif // DEBUG_PARAM_SETUP

            // CheckNoiseGrowth with current kappa and nu
            int omega = 1;//kappa /(2*(depth()+1));
            int tau = (sigma+omega-1)/omega;
            double err = 6.0;               // number of std devs; allows td times growth of noise reducing failure prob.
            double adds = 1.0;              // max number of additions per level before multiplication
            double noise[CIRCUIT_DEPTH];
            noise[0] = beta_;

            ComputeNoise(noise[0], nu, kappa, tau, omega, adds, err);
            for(int i=1; i<depth_; i++)
            {
                sigma -= kappa;
                tau = (sigma+omega-1)/omega;
                noise[i] = noise[i-1];
                ComputeNoise(noise[i], nu, kappa, tau, omega, adds, err);
            }

            #ifdef DEBUG_PARAM_SETUP
                cout << "Noise size after each level: ";
            #endif // DEBUG_PARAM_SETUP
            for(int i=0; i<depth_; i++)
            {
                noise[i] = log2(noise[i]);

                #ifdef DEBUG_PARAM_SETUP
                    cout << (int)noise[i] << " ";
                #endif // DEBUG_PARAM_SETUP
            }
            #ifdef DEBUG_PARAM_SETUP
                cout << endl;
            #endif // PARAM_SETUP

            if(noise[depth_-1] - noise[0] < noise[0])
                found = true;
            else
                kappa++;

            sigma = kappa * (depth()+1);
        }

    #ifdef DEBUG_PARAM_SETUP
        cout << endl << "AUTO PARAM SETUP END" << endl<< endl;
    #endif // DEBUG_PARAM_SETUP

        set_kappa(kappa);
        set_lambda(kappa);
        set_sigma(sigma);
        #ifdef XN_1
            set_nu(nu);
        #elif defined CYCLOTOMIC
            int mu = GenPrime_long(log2(nu));
            set_mu(mu);
            set_nu(EulerToient(mu));
        #endif // XN_1
        set_omega(kappa /(2*(depth()+1)));
        set_tau((sigma+omega_-1)/omega_);
        CheckSecurity();

    #endif // SECURE
}

void Params::ComputeNoise(double &noise, double nu, double kappa, double tau, double omega, double adds, double err)
{
    noise = err*( (adds*SqrRoot(nu)*noise*noise + pi_*nu*beta_*pow(2, omega)*((pi_+1)*beta_+1)*tau)/(pow(2,kappa)) + SqrRoot(nu)*pi_*(pi_*beta_+1));
}
//void Params::ComputeNoise(double n, double q, double k, double b, double td, double adds, double w)
//{
//    double result = td*(adds*SqrRoot(n)*pow(b,2) + pi_*n*b*w*((pi_+1)*b+1)*q )
    // noisea[0]=B0fa(nn,qq, BB,KK,tdd)
    // B0fa(n,q,K,B,mu)=mu*( (sqrt(n)*aa*B^2 + SM*n*B*w*((SM+1)*B+1)*lq)*K + sqrt(n)*SM*(SM*B+1) )
//}

double Params::ComputeHermite()
{
    double temp = ((double)sigma_-4)/(4*nu_);
    set_hermite(pow(2.0, temp));
    return hermite_;
}

bool Params::CheckSecurity()
{
    bool result = (ComputeHermite() < delta_);
#ifdef DEBUG_INFO
    cout << "Setup Secure : ";
    if(result) cout << "yes" << endl;
    else cout << "no" << endl;
#endif // DEBUG_INFO
    return result;
}

bool Params::CheckNoiseGrowth()
{

}

// SETs
void Params::set_setup_type(SetupType type)
{
    setup_type_ = type;
#ifdef DEBUG_INFO
    cout << "Setup Type : ";
    if(type == automatic) cout << "automatic" << endl;
    else cout << "manual" << endl;
#endif // DEBUG_INFO
}
void Params::set_hom_type(HomType type)
{
    hom_type_ = type;

#ifdef DEBUG_INFO
    cout << "Homomorphic Type : ";
    if(type == fhe) cout << "fully" << endl;
    else if(type == lhe) cout << "leveled" << endl;
    else cout << "somewhat" << endl;
#endif // DEBUG_INFO
}
void Params::set_ring_type(RingType type)
{
    ring_type_ = type;
#ifdef DEBUG_INFO
    cout << "Ring Type : ";
    if(type == cyclotomic) cout << "cyclotomic" << endl;
    else cout << "xn_1" << endl;
#endif // DEBUG_INFO
}
void Params::set_batch_flag(Flag flag)
{
    batch_flag_ = flag;
#ifdef DEBUG_INFO
    cout << "Batching : ";
    if(flag == on) cout << "on" << endl;
    else cout << "off" << endl;
#endif // DEBUG_INFO
}
void Params::set_fft_flag(Flag flag)
{
    fft_flag_ = flag;
#ifdef DEBUG_INFO
    cout << "FFT Mult : ";
    if(flag == on) cout << "on" << endl;
    else cout << "off" << endl;
#endif // DEBUG_INFO
}
void Params::set_reduc_type(ReductionType type)
{
    reduc_type_ = type;
#ifdef DEBUG_INFO
    cout << "Reduction Type : ";
    if(type == ntl) cout << "ntl" << endl;
    else if(type == fast) cout << "fast" << endl;
    else cout << "barrett" << endl;
#endif // DEBUG_INFO
}
void Params::set_pi(int val)
{
    pi_ = val;
#ifdef DEBUG_INFO
    cout << "Message Space p : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_mu(int val)
{
    mu_ = val;
#ifdef DEBUG_INFO
    cout << "Cyclotomic degree m : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_nu(int val)
{
    nu_ = val;
#ifdef DEBUG_INFO
    cout << "Ring poly degree n : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_omega(int val)
{
    omega_ = val;
#ifdef DEBUG_INFO
    cout << "Relin block bit-size w: " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_sigma(int val)
{
    sigma_ = val;
#ifdef DEBUG_INFO
    cout << "Max q bit-size : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_kappa(int val)
{
    kappa_ = val;
#ifdef DEBUG_INFO
    cout << "Cutting Size k : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_lambda(int val)
{
    lambda_ = val;
#ifdef DEBUG_INFO
    cout << "Min q bit-size : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_beta(int val)
{
    beta_ = val;
#ifdef DEBUG_INFO
    cout << "Initial noise bound B : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_depth(int val)
{
    depth_ = val;
#ifdef DEBUG_INFO
    cout << "Circuit Depth : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_tau(int val)
{
    tau_ = val;
#ifdef DEBUG_INFO
    cout << "Relin block/Eval Key count : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_rho(int val)
{
    rho_ = val;
#ifdef DEBUG_INFO
    cout << "Factor/Message slot count : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_theta(int val)
{
    theta_ = val;
#ifdef DEBUG_INFO
    cout << "Factor degree d : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_delta(double val)
{
    delta_ = val;
#ifdef DEBUG_INFO
    cout << "Hermite limit : " << val << endl;
#endif // DEBUG_INFO
}
void Params::set_hermite(double val)
{
    hermite_ = val;
#ifdef DEBUG_INFO
    cout << "Hermite factor : " << val << endl;
#endif // DEBUG_INFO
}
