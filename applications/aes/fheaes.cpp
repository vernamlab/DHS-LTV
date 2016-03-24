#include "fheaes.h"

FheAES::FheAES(){
	qIndex=0;
}

void FheAES::LTVKeySetUp(){

	ZZ q 	= to_ZZ("2");
	ZZ B 	= to_ZZ("1");
	GlobalParam gp;
	Set(gp, q, B, Modulus_M);

	n = new ltv(Modulus_M, Num_Primes, Max_Prime, 1, 16);
	n->Func_SetModulus();
	n->Func_ModulusFindRing(Num_Primes, Max_Prime, Dif_Prime, *n->Pointer_PolyModulus());
	n->Func_ComputeKeysRingRelinFFT(Num_Primes, 1);
	num_mul=0;
	num_relin=0;
}

void FheAES::SetMessageBits(unsigned char *mess){
	for(int i=0; i<16; i++)
		for(int j=0; j<8; j++)
			Bits[i*8+j] = (mess[i]>>j)%2;
}

void FheAES::EncryptMessage(){
	for(int i=0; i<128; i++)
		Bits[i] = n->Prim_Encrypt(Bits[i], 0);
}

void PrintChar(ZZX *tt, int num){
	long t[num*8];
	int mat[num];
	int r=0;
	for(int i=0; i<num*8; i++){
		t[i] = to_long(coeff(tt[i], 0));
	}

	for(int i=0; i<num; i++){
		r=0;
		for(int j=0; j<8; j++){
			r = (r<<1)+t[i*8+8-j-1];
		}
		mat[i] = r;
	}

	for(int i=0; i<num; i++){
		r = mat[(i%4)*4+i/4];
		if((r>>4) == 0)
			cout << 0;
		if(r==0)
			cout << 0;
		cout << uppercase << hex << r;
	}




	cout << endl;
}

void DecryptAndPrint(ltv *n, ZZX *Bits, int j){
	ZZX tt[128];
	for(int t=0; t<128; t++)
		tt[t] = n->Prim_Decrypt(Bits[t], j);
	PrintChar(tt, 16);
}

void FheAES::AESEncryption(){
	myTimer mt[10];
	AddRoundKey(0);
	for(int j=0; j<10; j++){
		cout << "Start Block\t" << j << endl;

		mt[j].Start();
///////////////////////////////////////////////////////////////
		for(int i=0; i<16; i++){
			AESSBox(&Bits[i*8], &Bits[i*8]);
			qIndex = qIndex-4;
		}
		qIndex = qIndex+4;

		if(j != 9)
			n->Prim_ReduceKeysLevel(4+j*4);

		FHEAESShiftRows();
		if(j != 9)
			FHEAESMixColumns();

		AddRoundKey(j+1);
///////////////////////////////////////////////////////////////

		mt[j].Stop();
		mt[j].ShowTime("Time:\t\t");
		cout << "Round My \t " << j << ":\t\t\t";
		DecryptAndPrint(n, Bits, 4*j+4);
		cout << "Round Org\t " << j << ":\t\t\t";
		ShowRoundResult(j);
	}

	double totaltime = 0.0;
	for(int i=0; i<10; i++)
		totaltime+=mt[i].GetTime();

	cout << "Total Time is:\t\t" << totaltime << endl;

}

void FheAES::AddRoundKey(int num){
	for(int i=0; i<128; i++)
		Bits[i] = Bits[i]+RoundKeys[num][i];
}

void FheAES::AESSBox(ZZX *bits_out, ZZX *bits_in){
	ZZX bit[8], bit2[8];
	Copy(bit, bits_in, 8);
	MatrixMultiply(bit, bit, InvX);
	GF28Inv(bit2, bit);
	MatrixMultiply(bit2, bit2, MX);
	Inv(bit2[0], bit2[0]);
	Inv(bit2[1], bit2[1]);
	Inv(bit2[5], bit2[5]);
	Inv(bit2[6], bit2[6]);

	Copy(bits_out, bit2, 8);
}

void FheAES::FHEAESShiftRows(){

	ZZX tp[32];

	for(int i = 0; i < 4; i ++)
		for(int j = 0; j < 8; j ++)
			tp[i*8 + j] = Bits[32+i*8+j];

	for(int j = 0; j < 8; j ++)
		Bits[32+ 24 + j] = tp[j];

	for(int j = 0; j < 8; j ++)
		Bits[32+ 0 + j] = tp[8+j];

	for(int j = 0; j < 8; j ++)
		Bits[32+ 8 + j] = tp[16+j];

	for(int j = 0; j < 8; j ++)
		Bits[32+ 16 + j] = tp[24+j];




	for(int i = 0; i < 4; i ++)
		for(int j = 0; j < 8; j ++)
			tp[i*8 + j] = Bits[64+i*8+j];

	for(int j = 0; j < 8; j ++)
		Bits[64+ 16 + j] = tp[j];

	for(int j = 0; j < 8; j ++)
		Bits[64+ 24 + j] = tp[8+j];

	for(int j = 0; j < 8; j ++)
		Bits[64+ 0 + j] = tp[16+j];

	for(int j = 0; j < 8; j ++)
		Bits[64+ 8 + j] = tp[24+j];




	for(int i = 0; i < 4; i ++)
		for(int j = 0; j < 8; j ++)
			tp[i*8 + j] = Bits[96+i*8+j];

	for(int j = 0; j < 8; j ++)
		Bits[96+ 8 + j] = tp[j];

	for(int j = 0; j < 8; j ++)
		Bits[96+ 16 + j] = tp[8+j];

	for(int j = 0; j < 8; j ++)
		Bits[96+ 24 + j] = tp[16+j];

	for(int j = 0; j < 8; j ++)
		Bits[96+ 0 + j] = tp[24+j];

}

void FheAES::FHEAESMixColumns(){

	ZZX tp[32];
/*
	for(int col = 0; col < 4; col ++)
	{
		for(int i = 0 ; i < 32 ; i ++)
			Copy(&tp[i], &Bits[col*32+i], 1);


		for(int row = 0 ; row < 4 ; row ++)
			GF2VectorMultiply(&Bits[row*8+col*32], tp, MixCol[row]);
	}
*/


	for(int col = 0; col < 4; col ++)
	{
		for(int i = 0 ; i < 32 ; i ++)
			Copy(&tp[i], &Bits[i/8*32+i%8 + col*8], 1);


		for(int row = 0 ; row < 4 ; row ++)
			GF2VectorMultiply(Bits + col*8 + 32*row, tp, MixCol[row]);
	}

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FheAES::Copy(ZZX *out, ZZX *in, int num){
	for(int i=0; i<num; i++)
		out[i] = in[i];
}

void FheAES::VectorMultiply(ZZX &out, ZZX* in, const unsigned char* a){
	out = 0;
	for(int i = 0; i < 8; i ++)
		if(a[i] == 1)
			out = out + in[7-i];//output.Add(output,output,(*b[7-i]));
}

void FheAES::MatrixMultiply(ZZX *out, ZZX *in, const unsigned char a[][8]){

	ZZX tp[8];
	for(int i = 0 ; i < 8 ; i ++)
		VectorMultiply(tp[i], in, a[7-i]);

	for(int i = 0 ; i < 8 ; i ++)
		out[i] = tp[i];
}

void FheAES::RelinBits(ZZX *bits, int num){
	for(int i=0; i<num; i++)
		n->Arith_PolyCoeffReduce(bits[i], bits[i], qIndex);

	for(int i=0; i<num; i++)
		bits[i]	= n->Prim_RelinRingFFT(bits[i], qIndex);

	for(int i=0; i<num; i++)
		n->Arith_PolyCoeffReduce(bits[i], bits[i], qIndex);

	num_relin+=num;
}

void FheAES::ModSwitchBits(ZZX *bits, int num){

	for(int i=0; i<num; i++){
		bits[i] = n->Prim_ModSwitch(bits[i], qIndex);
		n->Arith_PolyCoeffReduce(bits[i], bits[i], qIndex+1);
	}
}

void FheAES::GF28Inv(ZZX *bits_out, ZZX *bits_in){
	ZZX tp[20];

	ZZX tp2[8];
	for(int i=0; i<8; i++)
		tp2[i] = bits_in[i];


////////////////////////////////////////////////
	GF2Add(tp, bits_in, bits_in+4, 4);
	GF24Mul(tp+4,bits_in,bits_in+4);   //Mult

	RelinBits(tp+4, 4);
	ModSwitchBits(tp, 8);
	ModSwitchBits(tp2, 8);
	qIndex++;

////////////////////////////////////////////////

	GF24NSqr(tp+8,tp);
	GF2Add(tp+12,tp+4,tp+8,4);
	GF24Inv(tp+16,tp+12); //Mult

	RelinBits(tp+16, 4);
	ModSwitchBits(tp+16, 4);
	ModSwitchBits(tp2, 8);
	qIndex++;

////////////////////////////////////////////////

//	GF24Mul(bits_out,  bits_in+4,tp+16); //Mult
	GF24Mul(bits_out, tp2+4, tp+16); //Mult
	RelinBits(bits_out, 4);
	ModSwitchBits(bits_out, 4);

////////////////////////////////////////////////

//	GF24Mul(bits_out+4,bits_in,  tp+16); //Mult
	GF24Mul(bits_out+4, tp2,  tp+16); //Mult
	RelinBits(bits_out+4, 4);
	ModSwitchBits(bits_out+4, 4);
	qIndex++;

}

void FheAES::GF24Mul(ZZX* bits_out, ZZX* bits_in0, ZZX* bits_in1){
	ZZX tp[12];

	GF2Add(tp, bits_in0, bits_in0+2, 2);
	GF2Add(tp+2, bits_in1, bits_in1+2, 2);

	GF22Mul(tp+4, bits_in0+2, bits_in1+2);
	GF22Mul(tp+6, tp, tp+2);
	GF22Mul(tp+8, bits_in0, bits_in1);

	GF22N(tp+10, tp+6);

	GF2Add(bits_out, tp+10, tp+8, 2);
	GF2Add(bits_out+2, tp+10, tp+4, 2);

}

void FheAES::GF24Inv(ZZX* out, ZZX* in){
	ZZX tp[12];

	GF2Add(tp,   in,in+2,2);
	GF22Mul(tp+2,in,in+2);

		RelinBits(tp+2, 2);
		ModSwitchBits(tp+2, 2);
		qIndex++;

	GF22Sqr(tp+4,tp);
	GF22N(tp+6,tp+4);

	GF2Add(tp+8,tp+2,tp+6,2);

	GF22Inv(tp+10,tp+8);

	GF22Mul(out,in+2,tp+10);
	GF22Mul(out+2,in,tp+10);

}

void FheAES::GF22Inv(ZZX* out, ZZX* in){
	GF22Sqr(out,in);
}

void FheAES::GF22N(ZZX* bits_out, ZZX* bits_in)
{
	Add(bits_out[0], bits_in[0], bits_in[1]);
	bits_out[1] = bits_in[0];
}

void FheAES::GF22Sqr(ZZX* bits_out,ZZX* bits_in){
	ZZX tp[2];

	tp[0] = bits_in[0];
	tp[1] = bits_in[1];

	bits_out[0] = tp[1];
	bits_out[1] = tp[0];
}

void FheAES::GF22Mul(ZZX* bits_out, ZZX* bits_in0, ZZX *bits_in1){
	ZZX tp[5];

	Add(tp[0], bits_in0[0], bits_in0[1]);
	Add(tp[1], bits_in1[0], bits_in1[1]);

	Mul(tp[2], bits_in0[1], bits_in1[1]);
	Mul(tp[3], tp[0], tp[1]);
	Mul(tp[4], bits_in0[0], bits_in1[0]);

	Add(bits_out[0], tp[3], tp[4]);
	Add(bits_out[1], tp[3], tp[2]);

}

void FheAES::GF24NSqr(ZZX* bits_out, ZZX* bits_in){
	ZZX tp[4];

	GF2Add(tp, bits_in, bits_in+2, 2);
	GF22N(tp+2, bits_in);

	GF22Sqr(bits_out+2, tp);
	GF22Sqr(bits_out, tp+2);
}

void FheAES::GF2Add(ZZX *bits_out, ZZX *bits_in0, ZZX *bits_in1, int num){
	for(int i=0; i<num; i++)
		GF2Add(&bits_out[i], &bits_in0[i], &bits_in1[i]);
}

void FheAES::GF2Add(ZZX *bits_out, ZZX *bits_in0, ZZX *bits_in1){
	*bits_out = *bits_in0 + *bits_in1;
}

void FheAES::Add(ZZX &bits_out, ZZX &bits_in0, ZZX &bits_in1){
	bits_out = bits_in0 + bits_in1;
}

void FheAES::Mul(ZZX &bits_out, ZZX &bits_in0, ZZX &bits_in1){
	bits_out = bits_in0*bits_in1;
	num_mul++;
}

void FheAES::Inv(ZZX &bits_out, ZZX &bits_in){
	bits_out = bits_in+1;
}


void Print(ZZX *c, int num){
	for(int i=0; i<num; i++)
		if(c[i] == 0)
			cout << 0;
		else
			cout << 1;
		cout << endl;
}

void FheAES::GF2VectorMultiply(ZZX* out, ZZX* a , const unsigned char* b){
	ZZX tp[8];

	GF2MultiplyConst(tp, a, b[0]);
	for(int i = 0 ; i < 8 ; i ++)
		Copy(&out[i], &tp[i], 1);

	GF2MultiplyConst(tp, a +8, b[1]);
	GF2Add(out, out, tp, 8);
	GF2MultiplyConst(tp, a +16, b[2]);
	GF2Add(out, out, tp, 8);
	GF2MultiplyConst(tp, a +24, b[3]);
	GF2Add(out, out, tp, 8);

}

void FheAES::GF2MultiplyConst(ZZX* out, ZZX* a , const unsigned char b){
	if(b == 1){
		for (int i= 0; i < 8; i ++)
			Copy(&out[i], &a[i], 1);
	}
	if(b == 2){
		for (int i= 1; i < 8; i ++)
			Copy(&out[i], &a[i-1], 1);


		Copy(&out[0], &a[7], 1);
		Add(out[1], out[1], a[7]);
		Add(out[3], out[3], a[7]);
		Add(out[4], out[4], a[7]);
	}
	if(b == 3){
		for (int i= 0; i < 8; i ++)
			Copy(&out[i], &a[i], 1);


		for (int i= 1; i < 8; i ++)
			Add(out[i], out[i], a[i-1]);

		Add(out[0], out[0], a[7]);
		Add(out[1], out[1], a[7]);
		Add(out[3], out[3], a[7]);
		Add(out[4], out[4], a[7]);
	}
}

void FheAES::SetKeys(){

	for(int s=0; s<11; s++)
		Char2ZZX(key_rounds[s], RoundKeys[s]);

	for(int s=0; s<11; s++)
		for(int i=0; i<128; i++)
			RoundKeys[s][i] = n->Prim_Encrypt(RoundKeys[s][i], 0);

	for(int s=0; s<11; s++)
			for(int i=0; i<128; i++)
				n->Arith_CoeffReduce(RoundKeys[s][i], RoundKeys[s][i], 4*s);
}

void FheAES::Char2ZZX(const unsigned char *c, ZZX *k){
	for(int i=0; i<16; i++)
		for(int j=0; j<8; j++)
			k[i*8+j] = (c[i]>>j)%2;
}

void PrintResChar(const unsigned char *c){
	int r=0;
	for(int i=0; i<16; i++){
		r = (int)c[i];
		if((r>>4)==0)
			cout << 0;
		if(r==0)
			cout << 0;
		cout << uppercase << hex << r;
	}
	cout << endl;
}

void FheAES::ShowRoundResult(int r){
	PrintResChar(result_rounds[r]);
}

void Copy(ZZX &out, ZZX &in){
	out = in;
}


///////////////////////////////////////////////////////////////////////////////























