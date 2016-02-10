#ifndef DEF_H_
#define DEF_H_

//#define MANUAL_SETUP
#define AUTO_SETUP

#define CYCLOTOMIC
//#define XN_1

#define SECURE

//#define FHE
#define LHE
//#define SHE

#define FAST_REDUCTION
//#define BARRETT_REDUCTION
//#define NTL_REDUCTION

//#define FFT
//#define BATCH

#define DEBUG_INFO
#define DEBUG_PARAM_SETUP

//#define TIMING_INFO                               // General setup timings
//#define TIMING_DETAIL                             // Inner method timings


// APPLICATION PARAMS
#define PI                  1                                       // bit size of p : Message Space Z_p. p = 2^{pi}.
#define CIRCUIT_DEPTH       3                                       // d : Circuit Depth
#define BETA                1                                       // bit size of B : Initial Noise Bound, i.e. sample s,g,e,f' in CHI are all bounded by 2^Beta.


// NOISE PARAMS
//#ifdef AUTO_SETUP                                                   // Parameter Selection is on auto-pilot.

#ifdef MANUAL_SETUP
    #define KAPPA               24                                      // K : Noise Cutting Bit Size is Kappa
    #ifdef SHE
        #define LAMBDA          30                                      // L : Smallest (last level) q bit size for the ring R = Z_{q_{d}}/F(x) is Lambda.
    #else
        #define LAMBDA          KAPPA                                   // In case of leveled hom., the last step q is the same size as cutting size
    #endif // SHE

    #define SIGMA               (LAMBDA+(KAPPA*CIRCUIT_DEPTH))          // Do not change this. Max q size.

    #ifdef LHE
        #define OMEGA          (KAPPA/(2*(CIRCUIT_DEPTH+1)))            // Relinearization block size; i.e. Ctext will be seperated into OMEGA_ bits of chunks
                                                                        // If you want, you can set this value to any integer {1, 2, ...}
        #define TAU             ((SIGMA+OMEGA-1)/OMEGA)                 // Do not change this even if you modify omega.
    #endif // LHE


    // SECURITY PARAMS
    #ifdef CYCLOTOMIC
        #define MU              2049                                    // M : Cyclotomic degree of the ring modulus F(x); i.e. deg(F) = N = Phi(M)
    #else
        #define NU              257                                    // N: Degree of the ring modulus F(x) = x^n - 1
    #endif // CYCLOTOMIC
#endif // MANUAL_SETUP


// SECURITY PARAMS
#ifdef SECURE
    #define DELTA           1.0066                                  // Hermite Factor will be smaller than this value
                                                                    // 1.0066 achieves 80-bit security
                                                                    // Hermite Factor = pow(2, (log(q_{0})-4)/(4*N)) < 1.0066
#endif // SECURE




#endif // DEF_H_
