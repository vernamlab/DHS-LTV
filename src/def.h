#ifndef DEF_H_
#define DEF_H_

//////FOR PHI OF M
#define M_SIZE 51

#if M_SIZE == 51			// Toy setting n = 32
	#define Modulus_M	51
#elif M_SIZE == 8191		// n = 8190
	#define Modulus_M	8191
#elif M_SIZE == 21845		// n = 16384
	#define Modulus_M	21845
#elif M_SIZE == 32767		// n = 27000
	#define Modulus_M	32767
#elif M_SIZE == 65535		// n = 32768
	#define Modulus_M	65535
#endif

////FOR MODULUS
#define Num_Primes  41
#define Max_Prime	(1271)
#define Dif_Prime   31

#endif /* DEF_H_ */
