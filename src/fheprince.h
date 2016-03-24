#ifndef FHEPRINCE_H_
#define FHEPRINCE_H_

#include "def.h"
#include "ltv.h"
#include "fft_mult.h"
#include "general.h"

class FhePRINCE{
public:
	FhePRINCE();
	void LTVKeySetUp();
	void SetMessage(int *m);
	void SetKeys(int *k0, int *k1);
	void EncryptMessage();
	void EncryptKeys();
	void PRINCEEncryption();
	int num_mul, num_relin;
private:
	ltv *n;
//	myCRT c;

	int qIndex, levels;
	ZZX Bits[64], Key_0[64], Key_1[64];


	void PRINCE_enc(ZZX PLAINTEXT[64], ZZX K_0[64], ZZX K_1[64]);
	void addRoundKey(ZZX in[64], ZZX key[64]);
	void addRC(ZZX in[64], int round);
	void SBOX(ZZX in[64]);
	void _sbox(ZZX in[4]);
	void C_AND(ZZX &r, ZZX &a, ZZX &b);
	void MixColumn(ZZX in[64]);
	void M_p(ZZX in[64]);
	void ShiftRow(ZZX in[64]);
	void INV_SBOX(ZZX in[64]);
	void _inv_sbox(ZZX in[4]);
	void inv_MixColumn(ZZX in[64]);
	void inv_ShiftRow(ZZX in[64]);
	void KeyExpansion(ZZX key[64]);

};



#endif /* FHEPRINCE_H_ */
