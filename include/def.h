#ifndef DEF_H_
#define DEF_H_

// HOMOMORPHIC TYPE
//#define SHE
#define FHE
#define AUTO_PARAM_SETUP

// APPLICATION PARAMS
#define P 2                                        // Message Space Z_p
#define DEPTH 1                                    // Circuit Depth

// NOISE PARAMS
#define KAPPA 20                                   // Noise Cutting Bit Size
#ifdef SHE
    #define LAMBDA 30                              // Smallest (last level) q bit size for the ring R = Z_{q_{d}}/F(x)
#else
    #define LAMBDA KAPPA                          // In case leveled FHE it is the same as Kappa
#endif // SHE
#define SIGMA  (LAMBDA + (KAPPA*DEPTH))         // Largest (first level) q bit size for the ring R = Z_{q_{0}}/F(x)
                                                    // SIGMA_ = log_2(q_{0})
#ifdef FHE
    #define OMEGA (KAPPA/(2*(DEPTH+1)))          // Relinearization block size; i.e. Ctext will be seperated into OMEGA_ bits of chunks
    #define TAU ((SIGMA+OMEGA-1)/OMEGA)         // Eval Keys will be : E(f*2^{OMEGA_*i}) where i=0,...,TAU_-1
#endif // FHE

// SECURITY PARAMS
#define N 257                                      // Degree of the ring poly : R = Z_{q_{i}}/F(x), for each i=0,...,d
                                                    // N = deg(F)
#ifdef AUTO_PARAM_SETUP                             // Parameter Selection is on auto-pilot.
    #define DELTA 1.0066                           // Hermite Factor will be smaller than this value
                                                    // 1.0066 achieves 80-bit security
                                                    // Hermite Factor = pow(2, (log(q_{0})-4)/(4*N)) < 1.0066
#endif // AUTO_PARAM_SETUP


// POLYNOMIAL SETUP
#define RING_TYPE       cyclotomic                      // {cyclotomic, xn_1}
#define REDUCTION_TYPE  fast                            // {fast, barrett, ntl}
#define FFT             off                             // {on, off}
#define BATCH           off                             // {on, off}



#endif // DEF_H_
