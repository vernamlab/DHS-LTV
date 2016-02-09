#ifndef DEF_H_
#define DEF_H_



// POLYNOMIAL SETUP
#define HOMO_TYPE           fhe                             // {fhe, she}
#define RING_TYPE           cyclotomic                      // {cyclotomic, xn_1}
#define REDUCTION_TYPE      fast                            // {fast, barrett, ntl}
#define FLAG_FFT            off                             // {on, off}
#define FLAG_BATCH          off                             // {on, off}
#define PARAM_SETUP         manual                          // {auto, manual}



// APPLICATION PARAMS
#define MSG_DOMAIN          2                               // p : Message Space Z_p
#define CIRCUIT_DEPTH       1                               // d : Circuit Depth



// NOISE PARAMS
#define BETA                1                               // B : Initial Noise Bound, i.e. sample s,g,e,f' in CHI are all bounded by Beta.
#define KAPPA               20                              // K : Noise Cutting Bit Size is Kappa
#if HOMO_TYPE == she
    #define LAMBDA          30                              // L : Smallest (last level) q bit size for the ring R = Z_{q_{d}}/F(x) is Lambda.
#else
    #define LAMBDA          KAPPA                           // L = K : In case leveled FHE it is the same as Kappa
#endif // HOMO_TYPE
#define SIGMA               (LAMBDA + (KAPPA*DEPTH))        // Largest (first level) q bit size for the ring R = Z_{q_{0}}/F(x)
                                                            // SIGMA_ = log_2(q_{0})
#if HOMO_TYPE == fhe
    #define OMEGA           (KAPPA/(2*(DEPTH+1)))           // Relinearization block size; i.e. Ctext will be seperated into OMEGA_ bits of chunks
    #define TAU             ((SIGMA+OMEGA-1)/OMEGA)         // Eval Keys will be : E(f*2^{OMEGA_*i}) where i=0,...,TAU_-1
#endif // HOMO_TYPE



// SECURITY PARAMS
#if PARAM_SETUP == auto                                     // Parameter Selection is on auto-pilot.
    #define DELTA           1.0066                          // Hermite Factor will be smaller than this value
#endif // PARAM_SETUP                                       // 1.0066 achieves 80-bit security
                                                            // Hermite Factor = pow(2, (log(q_{0})-4)/(4*N)) < 1.0066
#if RING_TYPE == cyclotomic
    #define MU              257                             // M : Cyclotomic degree of the ring modulus F(x); i.e. deg(F) = N = Phi(M)
#else
    #define NU              257                             // N: Degree of the ring modulus F(x) = x^n - 1
#endif // RING_TYPE



#endif // DEF_H_
