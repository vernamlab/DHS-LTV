#ifndef DEF_H_
#define DEF_H_

//Randomness
#if Random == 1
	#define PureRandom
#else
	#define PseudoRandom
#endif

typedef enum SortingAlgorithm
{
    Bubble = 1,
    Batcher = 2,
    Direct = 3,
    Greedy = 4
} SortingAlgorithm;

#define Sort_Alg        1
#define Sort_Size       4
#define Bit_Size        8


//////////////////////////////////////

#if Bit_Size == 8
    #define LT_Depth    4
#elif Bit_Size == 32
    #define LT_Depth    6
#endif // Bit_Size

#define CS_Depth        (LT_Depth+1)

#define EQ_Depth        3

///////////////////////////////////////////////////////
#if Sort_Alg == 1
    #define Depth               (CS_Depth*Sort_Size)

#elif Sort_Alg == 2
    #if Sort_Size == 4
        #define Other_Depth     (3)
    #elif Sort_Size == 8
        #define Other_Depth     (6)
    #elif Sort_Size == 16
        #define Other_Depth     (10)
    #elif Sort_Size == 32
        #define Other_Depth     (15)
    #elif Sort_Size == 64
        #define Other_Depth     (36)
    #endif // Sort_Size

    #define Depth               (CS_Depth*Other_Depth)

#elif Sort_Alg == 3
    #if Sort_Size == 4
        #define Other_Depth     (1)
    #elif Sort_Size == 8
        #define Other_Depth     (2)
    #elif Sort_Size == 16
        #define Other_Depth     (3)
    #elif Sort_Size == 32
        #define Other_Depth     (4)
    #elif Sort_Size == 64
        #define Other_Depth     (5)
    #endif // Sort_Size

    #define Depth               (LT_Depth+Other_Depth+EQ_Depth+1)

#elif Sort_Alg == 4
    #if Sort_Size == 4
        #define Other_Depth     (2)
    #elif Sort_Size == 8
        #define Other_Depth     (3)
    #elif Sort_Size == 16
        #define Other_Depth     (4)
    #elif Sort_Size == 32
        #define Other_Depth     (5)
    #elif Sort_Size == 64
        #define Other_Depth     (6)
    #endif // Sort_Size

    #define Depth               (LT_Depth+Other_Depth+1)

#endif // Sort_Alg
///////////////////////////////////////////////////////

#define Num_Primes              (Depth+1)
#define Dif_Prime               20                              // Increase, if the noise growth is large
#define Max_Prime	            (Num_Primes*Dif_Prime)
#define Modulus_M               15                              // Poly modulus, N = Phi(Modulus_M)
                                                                // Message Slot Size = N/FactorDeg
                                                                // Ex; 630 = Phi(8191)/13 = 8190/13
                                                                // For test purposes, it is better to keep degree m as low as possible.
                                                                // For security purposes, it is better to keep degree m as high as possible.

// Popular M values:
// Modulus_M 15
// Modulus_M 255
// Modulus_M 8191
// Modulus_M 21845
// Modulus_M 32767



#endif /* DEF_H_ */
