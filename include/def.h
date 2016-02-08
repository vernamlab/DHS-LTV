#ifndef DEF_H_
#define DEF_H_

// HOMOMORPHIC TYPE
#define SHE_
//#define FHE_

// APPLICATION PARAMS
#define P_ 2                            // Message Space Z_p
#define DEPTH_ 1                        // Circuit Depth

// NOISE PARAMS
#define KAPPA_ 20                                   // Noise Cutting Bit Size
#ifdef SHE_
    #define LAMBDA_ 30                              // Smallest (last level) q bit size for the ring R = Z_{q_{d}}/F(x)
#else
    #define LAMBDA_ KAPPA_                          // In case leveled FHE it is the same as Kappa
#endif // SHE_
#define MAXq_  (LAMBDA_ + (KAPPA_*DEPTH_))          // Largest (first level) q bit size for the ring R = Z_{q_{0}}/F(x)
                                                    // MAXq_ = log_2(q_{0})

// SECURITY PARAMS
#define N_ 257                          // Degree of the ring poly : R = Z_{q_{i}}/F(x), for each i=0,...,d
                                        // N = deg(F)
#define DELTA_ 1.0066                   // Hermite Factor needs to be smaller than this value
                                        // 1.0066 achieves 80-bit security
                                        // Hermite Factor = pow(2, (log(q_{0})-4)/(4*N)) < 1.0066


#endif // DEF_H_
