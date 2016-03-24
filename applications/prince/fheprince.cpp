#include "fheprince.h"

FhePRINCE::FhePRINCE(){
	qIndex = 0;
	levels = 0;
	num_mul = 0;
	num_relin = 0;
}

void FhePRINCE::LTVKeySetUp(){

	ZZ q 	= to_ZZ("2");
	ZZ B 	= to_ZZ("1");
	GlobalParam gp;
	Set(gp, q, B, Modulus_M);

	n = new ltv(Modulus_M, Num_Primes, Max_Prime, 1, 16);
	n->Func_SetModulus();
	n->Func_ModulusFindRing(Num_Primes, Max_Prime, Dif_Prime, *n->Pointer_PolyModulus());
	n->Func_ComputeKeysRingRelinFFT(Num_Primes, 1);
}

void FhePRINCE::SetMessage(int *m){
	for(int i=0; i<64; i++)
		Bits[i] = m[i];
}

void FhePRINCE::SetKeys(int *k0, int *k1){
	for(int i=0; i<64; i++){
		Key_0[i] = k0[i];
		Key_1[i] = k1[i];
	}
}

void FhePRINCE::EncryptMessage(){
	for(int i=0; i<64; i++)
		Bits[i] = n->Prim_Encrypt(Bits[i], 0);
}

void FhePRINCE::EncryptKeys(){
	for(int i=0; i<64; i++){
		Key_0[i] = n->Prim_Encrypt(Key_0[i], 0);
		Key_1[i] = n->Prim_Encrypt(Key_1[i], 0);
	}
}

void FhePRINCE::PRINCEEncryption(){

	myTimer t;
	t.Start();
	PRINCE_enc(Bits, Key_0, Key_1);
	t.Stop();
	t.ShowTime("PRINCE Time:\t");
	ZZX res[64];
	for(int i=0; i<64; i++)
		res[i] = n->Prim_Decrypt(Bits[i], Num_Primes-1);

	cout << endl << endl;
	for(int i=0; i<64; i++)
		cout << coeff(res[i], 0);

	cout << endl;
	cout << "1001111110110101000110010011010111111100001111011111010100100100" << endl;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int RC[768]=
    {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,1,1,0,0,0,1,1,0,0,1,1,0,0,0,1,0,1,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,1,1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,0,0,1,0,0,
    1,0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,1,1,0,0,1,1,1,1,1,0,0,1,1,0,0,0,1,1,1,0,1,0,0,0,0,
    0,0,0,0,1,0,0,0,0,0,1,0,1,1,1,0,1,1,1,1,1,0,1,0,1,0,0,1,1,0,0,0,1,1,1,0,1,1,0,0,0,1,0,0,1,1,1,0,0,1,1,0,1,1,0,0,1,0,0,0,1,0,0,1,
    0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,1,1,0,0,0,1,1,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,1,0,0,1,1,0,1,1,1,0,1,1,1,
    1,0,1,1,1,1,1,0,0,1,0,1,0,1,0,0,0,1,1,0,0,1,1,0,1,1,0,0,1,1,1,1,0,0,1,1,0,1,0,0,1,1,1,0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,1,0,1,1,0,0,
    0,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0,0,1,0,0,1,1,1,1,0,1,1,1,1,0,0,0,1,1,1,1,1,1,0,1,1,0,0,1,0,1,0,1,0,1,0,1,1,1,0,0,1,0,1,1,0,0,0,1,
    1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,1,1,1,1,1,0,0,0,1,1,0,1,0,1,1,0,0,0,1,0,0,0,0,1,1,1,0,1,0,1,0,1,0,
    1,1,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,1,1,0,0,1,0,0,1,0,1,0,0,1,1,0,0,1,0,0,0,1,1,1,1,0,0,0,1,0,1,0,1,0,0,
    0,1,1,0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0,1,0,1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,1,1,0,1,
    1,1,0,1,0,0,1,1,1,0,1,1,0,1,0,1,1,0,1,0,0,0,1,1,1,0,0,1,1,0,0,1,1,1,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,1,1,1,0,0,1,1,0,0,1,
    1,1,0,0,0,0,0,0,1,0,1,0,1,1,0,0,0,0,1,0,1,0,0,1,1,0,1,1,0,1,1,1,1,1,0,0,1,0,0,1,0,1,1,1,1,1,0,0,0,1,0,1,0,0,0,0,1,1,0,1,1,1,0,1,
    };
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FhePRINCE::addRoundKey(ZZX in[64], ZZX key[64]){
    for (int i=0; i<64; i++)
        in[i] = in[i] + key[i];
}

void FhePRINCE::addRC(ZZX in[64], int round){
    for (int i=0;i<64;i++)
        in[i] = (in[i] + RC[(round*64)+i]);
}

void FhePRINCE::C_AND(ZZX &r, ZZX &a, ZZX &b){
    r = a*b;
    num_mul++;
}

void FhePRINCE::_sbox(ZZX in[4])
{
    ZZX A, B, C, D;
    ZZX AB,AC,AD,BC,BD,CD;
    ZZX ABD,ACD,BCD,ABC;

    A = in[0];
    B = in[1];
    C = in[2];
    D = in[3];

    C_AND(AB, A, B);
    C_AND(AC, A, C);
    C_AND(AD, A, D);
    C_AND(BC, B, C);
    C_AND(BD, B, D);
    C_AND(CD, C, D);

    n->Arith_PolyCoeffReduce(AB, AB, qIndex);
    n->Arith_PolyCoeffReduce(AC, AC, qIndex);
    n->Arith_PolyCoeffReduce(AD, AD, qIndex);
    n->Arith_PolyCoeffReduce(BC, BC, qIndex);
    n->Arith_PolyCoeffReduce(BD, BD, qIndex);
    n->Arith_PolyCoeffReduce(CD, CD, qIndex);


    AB	= n->Prim_RelinRingFFT(AB, qIndex);
    CD	= n->Prim_RelinRingFFT(CD, qIndex);
    num_relin+=2;
    n->Arith_PolyCoeffReduce(AB, AB, qIndex);
    n->Arith_PolyCoeffReduce(CD, CD, qIndex);


    AB = n->Prim_ModSwitch(AB, qIndex);
    AC = n->Prim_ModSwitch(AC, qIndex);
    AD = n->Prim_ModSwitch(AD, qIndex);
    BC = n->Prim_ModSwitch(BC, qIndex);
    BD = n->Prim_ModSwitch(BD, qIndex);
    CD = n->Prim_ModSwitch(CD, qIndex);
    A = n->Prim_ModSwitch(A, qIndex);
    B = n->Prim_ModSwitch(B, qIndex);
    C = n->Prim_ModSwitch(C, qIndex);
    D = n->Prim_ModSwitch(D, qIndex);

    qIndex++;



    C_AND(ABD, AB, D);
    C_AND(ACD, CD, A);
    C_AND(BCD, CD, B);
    C_AND(ABC, AB, C);
    //(A+C+AB+BC+ABD+ACD+BCD)'
    in[0] = (A + C + AB + BC + ABD + ACD + BCD + 1);
    //A+D+AC+AD+CD+ABC+ACD
    in[1] = (A + D + AC + AD + CD + ABC + ACD);
    //(AC+BC+BD+ABC+BCD)'
    in[2] = (AC + BC + BD + ABC + BCD + 1);
    //(A+B+AB+AD+BC+CD+BCD)'
    in[3] = (A + B + AB + AD + BC + CD + BCD + 1);

    n->Arith_PolyCoeffReduce(in[0], in[0], qIndex);
    n->Arith_PolyCoeffReduce(in[1], in[1], qIndex);
    n->Arith_PolyCoeffReduce(in[2], in[2], qIndex);
    n->Arith_PolyCoeffReduce(in[3], in[3], qIndex);
    in[0]	= n->Prim_RelinRingFFT(in[0], qIndex);
    in[1]	= n->Prim_RelinRingFFT(in[1], qIndex);
    in[2]	= n->Prim_RelinRingFFT(in[2], qIndex);
    in[3]	= n->Prim_RelinRingFFT(in[3], qIndex);
    num_relin+=4;
    n->Arith_PolyCoeffReduce(in[0], in[0], qIndex);
    n->Arith_PolyCoeffReduce(in[1], in[1], qIndex);
    n->Arith_PolyCoeffReduce(in[2], in[2], qIndex);
    n->Arith_PolyCoeffReduce(in[3], in[3], qIndex);

    in[0] = n->Prim_ModSwitch(in[0], qIndex);
    in[1] = n->Prim_ModSwitch(in[1], qIndex);
    in[2] = n->Prim_ModSwitch(in[2], qIndex);
    in[3] = n->Prim_ModSwitch(in[3], qIndex);

    qIndex++;
}

void FhePRINCE::SBOX(ZZX in[64]){

    for (int i=0;i<16;i++){
        ZZX A[4] = {in[4*i],in[4*i+1],in[4*i+2],in[4*i+3]};
        _sbox(A);
        in[4*i] = A[0];
        in[4*i+1] = A[1];
        in[4*i+2] = A[2];
        in[4*i+3] = A[3];
        qIndex = qIndex-2;
    }
    qIndex = qIndex + 2;
    levels = levels + 2;
    n->Prim_ReduceKeysLevel(levels);
//  n->ReduceKeysLevel(levels+1, 1);
}

void FhePRINCE::_inv_sbox(ZZX in[4])
{
	ZZX A, B, C, D;
	ZZX AB,AC,AD,BC,BD,CD;
	ZZX ABC,ABD,BCD,ACD;

    A = in[0];
    B = in[1];
    C = in[2];
    D = in[3];

    C_AND(AB, A, B);
    C_AND(AC, A, C);
    C_AND(AD, A, D);
    C_AND(BC, B, C);
    C_AND(BD, B, D);
    C_AND(CD, C, D);

    n->Arith_PolyCoeffReduce(AB, AB, qIndex);
    n->Arith_PolyCoeffReduce(AC, AC, qIndex);
    n->Arith_PolyCoeffReduce(AD, AD, qIndex);
    n->Arith_PolyCoeffReduce(BC, BC, qIndex);
    n->Arith_PolyCoeffReduce(BD, BD, qIndex);
    n->Arith_PolyCoeffReduce(CD, CD, qIndex);

    AB	= n->Prim_RelinRingFFT(AB, qIndex, qIndex%2);
    CD	= n->Prim_RelinRingFFT(CD, qIndex, qIndex%2);
    num_relin+=2;
    n->Arith_PolyCoeffReduce(AB, AB, qIndex);
    n->Arith_PolyCoeffReduce(CD, CD, qIndex);

    AB = n->Prim_ModSwitch(AB, qIndex);
    AC = n->Prim_ModSwitch(AC, qIndex);
    AD = n->Prim_ModSwitch(AD, qIndex);
    BC = n->Prim_ModSwitch(BC, qIndex);
    BD = n->Prim_ModSwitch(BD, qIndex);
    CD = n->Prim_ModSwitch(CD, qIndex);

    A = n->Prim_ModSwitch(A, qIndex);
    B = n->Prim_ModSwitch(B, qIndex);
    C = n->Prim_ModSwitch(C, qIndex);
    D = n->Prim_ModSwitch(D, qIndex);

    qIndex++;

    C_AND(ABC, AB, C);
    C_AND(ABD, AB, D);
    C_AND(BCD, CD, B);
    C_AND(ACD, CD, A);

    //(C+D+AB+BC+BD+CD+ABC+ABD+BCD)'
    in[0] = (C + D + AB + BC + BD + CD + ABC + ABD + BCD + 1);
    //B+D+AC+BC+BD+CD+ACD+BCD
    in[1] = (B + D + AC + BC + BD + CD + ACD + BCD);
    //(AB+AC+BC+BD+BCD)'
    in[2] = (AB + AC + BC + BD + BCD + 1);
    //(A+AB+BC+CD+ABD+ACD)'
    in[3] = (A + AB + BC + CD + ABD + ACD + 1);

    n->Arith_PolyCoeffReduce(in[0], in[0], qIndex);
    n->Arith_PolyCoeffReduce(in[1], in[1], qIndex);
    n->Arith_PolyCoeffReduce(in[2], in[2], qIndex);
    n->Arith_PolyCoeffReduce(in[3], in[3], qIndex);

    in[0]	= n->Prim_RelinRingFFT(in[0], qIndex);
    in[1]	= n->Prim_RelinRingFFT(in[1], qIndex);
    in[2]	= n->Prim_RelinRingFFT(in[2], qIndex);
    in[3]	= n->Prim_RelinRingFFT(in[3], qIndex);
    num_relin+=4;
    n->Arith_PolyCoeffReduce(in[0], in[0], qIndex);
    n->Arith_PolyCoeffReduce(in[1], in[1], qIndex);
    n->Arith_PolyCoeffReduce(in[2], in[2], qIndex);
    n->Arith_PolyCoeffReduce(in[3], in[3], qIndex);

    in[0] = n->Prim_ModSwitch(in[0], qIndex);
    in[1] = n->Prim_ModSwitch(in[1], qIndex);
    in[2] = n->Prim_ModSwitch(in[2], qIndex);
    in[3] = n->Prim_ModSwitch(in[3], qIndex);

    qIndex++;
}

void FhePRINCE::INV_SBOX(ZZX in[64]){
    for (int i=0;i<16;i++){
        ZZX A[4] = {in[4*i],in[4*i+1],in[4*i+2],in[4*i+3]};
        _inv_sbox(A);
        in[4*i] = A[0];
        in[4*i+1] = A[1];
        in[4*i+2] = A[2];
        in[4*i+3] = A[3];
        qIndex = qIndex-2;
    }
    qIndex = qIndex + 2;
    levels = levels + 2;
    n->Prim_ReduceKeysLevel(levels);
//    n->ReduceKeysLevel(levels+1, 1);
}

void FhePRINCE::M_p(ZZX in[64]){
    ZZX in1[64];

    in1[0] = (in[4]+in[8]+in[12]);
    in1[1] = (in[1]+in[9]+in[13]);
    in1[2] = (in[2]+in[6]+in[14]);
    in1[3] = (in[3]+in[7]+in[11]);
    in1[4] = (in[0]+in[4]+in[8]);
    in1[5] = (in[5]+in[9]+in[13]);
    in1[6] = (in[2]+in[10]+in[14]);
    in1[7] = (in[3]+in[7]+in[15]);
    in1[8] = (in[0]+in[4]+in[12]);
    in1[9] = (in[1]+in[5]+in[9]);
    in1[10] = (in[6]+in[10]+in[14]);
    in1[11] = (in[3]+in[11]+in[15]);
    in1[12] = (in[0]+in[8]+in[12]);
    in1[13] = (in[1]+in[5]+in[13]);
    in1[14] = (in[2]+in[6]+in[10]);
    in1[15] = (in[7]+in[11]+in[15]);

    in1[16] = (in[16]+in[20]+in[24]);
    in1[17] = (in[21]+in[25]+in[29]);
    in1[18] = (in[18]+in[26]+in[30]);
    in1[19] = (in[19]+in[23]+in[31]);
    in1[20] = (in[16]+in[20]+in[28]);
    in1[21] = (in[17]+in[21]+in[25]);
    in1[22] = (in[22]+in[26]+in[30]);
    in1[23] = (in[19]+in[27]+in[31]);
    in1[24] = (in[16]+in[24]+in[28]);
    in1[25] = (in[17]+in[21]+in[29]);
    in1[26] = (in[18]+in[22]+in[26]);
    in1[27] = (in[23]+in[27]+in[31]);
    in1[28] = (in[20]+in[24]+in[28]);
    in1[29] = (in[17]+in[25]+in[29]);
    in1[30] = (in[18]+in[22]+in[30]);
    in1[31] = (in[19]+in[23]+in[27]);

    in1[32] = (in[32]+in[36]+in[40]);
    in1[33] = (in[37]+in[41]+in[45]);
    in1[34] = (in[34]+in[42]+in[46]);
    in1[35] = (in[35]+in[39]+in[47]);
    in1[36] = (in[32]+in[36]+in[44]);
    in1[37] = (in[33]+in[37]+in[41]);
    in1[38] = (in[38]+in[42]+in[46]);
    in1[39] = (in[35]+in[43]+in[47]);
    in1[40] = (in[32]+in[40]+in[44]);
    in1[41] = (in[33]+in[37]+in[45]);
    in1[42] = (in[34]+in[38]+in[42]);
    in1[43] = (in[39]+in[43]+in[47]);
    in1[44] = (in[36]+in[40]+in[44]);
    in1[45] = (in[33]+in[41]+in[45]);
    in1[46] = (in[34]+in[38]+in[46]);
    in1[47] = (in[35]+in[39]+in[43]);

    in1[48] = (in[52]+in[56]+in[60]);
    in1[49] = (in[49]+in[57]+in[61]);
    in1[50] = (in[50]+in[54]+in[62]);
    in1[51] = (in[51]+in[55]+in[59]);
    in1[52] = (in[48]+in[52]+in[56]);
    in1[53] = (in[53]+in[57]+in[61]);
    in1[54] = (in[50]+in[58]+in[62]);
    in1[55] = (in[51]+in[55]+in[63]);
    in1[56] = (in[48]+in[52]+in[60]);
    in1[57] = (in[49]+in[53]+in[57]);
    in1[58] = (in[54]+in[58]+in[62]);
    in1[59] = (in[51]+in[59]+in[63]);
    in1[60] = (in[48]+in[56]+in[60]);
    in1[61] = (in[49]+in[53]+in[61]);
    in1[62] = (in[50]+in[54]+in[58]);
    in1[63] = (in[55]+in[59]+in[63]);

    for (int i=0;i<64;i++)
        in[i] = in1[i];
}

void FhePRINCE::ShiftRow(ZZX in[64]){
    int i = 4;
    ZZX temp[16] = {in[i],in[i+1],in[i+2],in[i+3],in[i+16],in[i+17],in[i+18],in[i+19],in[i+32],in[i+33],in[i+34],in[i+35],in[i+48],in[i+49],in[i+50],in[i+51]};
    in[i] = temp[4];
    in[i+1] = temp[5];
    in[i+2] = temp[6];
    in[i+3] = temp[7];
    in[i+16] = temp[8];
    in[i+17] = temp[9];
    in[i+18] = temp[10];
    in[i+19] = temp[11];
    in[i+32] = temp[12];
    in[i+33] = temp[13];
    in[i+34] = temp[14];
    in[i+35] = temp[15];
    in[i+48] = temp[0];
    in[i+49] = temp[1];
    in[i+50] = temp[2];
    in[i+51] = temp[3];

    i = 8;
    ZZX temp1[16] = {in[i],in[i+1],in[i+2],in[i+3],in[i+16],in[i+17],in[i+18],in[i+19],in[i+32],in[i+33],in[i+34],in[i+35],in[i+48],in[i+49],in[i+50],in[i+51]};
    in[i] = temp1[8];
    in[i+1] = temp1[9];
    in[i+2] = temp1[10];
    in[i+3] = temp1[11];
    in[i+16] = temp1[12];
    in[i+17] = temp1[13];
    in[i+18] = temp1[14];
    in[i+19] = temp1[15];
    in[i+32] = temp1[0];
    in[i+33] = temp1[1];
    in[i+34] = temp1[2];
    in[i+35] = temp1[3];
    in[i+48] = temp1[4];
    in[i+49] = temp1[5];
    in[i+50] = temp1[6];
    in[i+51] = temp1[7];

    i = 12;
    ZZX temp2[16] = {in[i],in[i+1],in[i+2],in[i+3],in[i+16],in[i+17],in[i+18],in[i+19],in[i+32],in[i+33],in[i+34],in[i+35],in[i+48],in[i+49],in[i+50],in[i+51]};
    in[i] = temp2[12];
    in[i+1] = temp2[13];
    in[i+2] = temp2[14];
    in[i+3] = temp2[15];
    in[i+16] = temp2[0];
    in[i+17] = temp2[1];
    in[i+18] = temp2[2];
    in[i+19] = temp2[3];
    in[i+32] = temp2[4];
    in[i+33] = temp2[5];
    in[i+34] = temp2[6];
    in[i+35] = temp2[7];
    in[i+48] = temp2[8];
    in[i+49] = temp2[9];
    in[i+50] = temp2[10];
    in[i+51] = temp2[11];
}

void FhePRINCE::MixColumn(ZZX in[64]){
    M_p(in);
    ShiftRow(in);
}

void FhePRINCE::inv_ShiftRow(ZZX in[64]){
    int i = 4;
    ZZX temp2[16] = {in[i],in[i+1],in[i+2],in[i+3],in[i+16],in[i+17],in[i+18],in[i+19],in[i+32],in[i+33],in[i+34],in[i+35],in[i+48],in[i+49],in[i+50],in[i+51]};
    in[i] = temp2[12];
    in[i+1] = temp2[13];
    in[i+2] = temp2[14];
    in[i+3] = temp2[15];
    in[i+16] = temp2[0];
    in[i+17] = temp2[1];
    in[i+18] = temp2[2];
    in[i+19] = temp2[3];
    in[i+32] = temp2[4];
    in[i+33] = temp2[5];
    in[i+34] = temp2[6];
    in[i+35] = temp2[7];
    in[i+48] = temp2[8];
    in[i+49] = temp2[9];
    in[i+50] = temp2[10];
    in[i+51] = temp2[11];

    i = 8;
    ZZX temp1[16] = {in[i],in[i+1],in[i+2],in[i+3],in[i+16],in[i+17],in[i+18],in[i+19],in[i+32],in[i+33],in[i+34],in[i+35],in[i+48],in[i+49],in[i+50],in[i+51]};
    in[i] = temp1[8];
    in[i+1] = temp1[9];
    in[i+2] = temp1[10];
    in[i+3] = temp1[11];
    in[i+16] = temp1[12];
    in[i+17] = temp1[13];
    in[i+18] = temp1[14];
    in[i+19] = temp1[15];
    in[i+32] = temp1[0];
    in[i+33] = temp1[1];
    in[i+34] = temp1[2];
    in[i+35] = temp1[3];
    in[i+48] = temp1[4];
    in[i+49] = temp1[5];
    in[i+50] = temp1[6];
    in[i+51] = temp1[7];

    i = 12;
    ZZX temp[16] = {in[i],in[i+1],in[i+2],in[i+3],in[i+16],in[i+17],in[i+18],in[i+19],in[i+32],in[i+33],in[i+34],in[i+35],in[i+48],in[i+49],in[i+50],in[i+51]};
    in[i] = temp[4];
    in[i+1] = temp[5];
    in[i+2] = temp[6];
    in[i+3] = temp[7];
    in[i+16] = temp[8];
    in[i+17] = temp[9];
    in[i+18] = temp[10];
    in[i+19] = temp[11];
    in[i+32] = temp[12];
    in[i+33] = temp[13];
    in[i+34] = temp[14];
    in[i+35] = temp[15];
    in[i+48] = temp[0];
    in[i+49] = temp[1];
    in[i+50] = temp[2];
    in[i+51] = temp[3];
}

void FhePRINCE::inv_MixColumn(ZZX in[64]){
    inv_ShiftRow(in);
    M_p(in);
}

void FhePRINCE::KeyExpansion(ZZX key[64]){
    ZZX temp[64];
    for (int i=0;i<64;i++)
        temp[i] = key[i];
    key[0] = temp[63];
    for (int i=0;i<63;i++)
        key[i+1] = temp[i];
    key[63] = (key[63]+temp[0]);
}

void FhePRINCE::PRINCE_enc(ZZX PLAINTEXT[64], ZZX K_0[64], ZZX K_1[64]){

    int ROUND = 0;
    addRoundKey(PLAINTEXT, K_0);
    addRoundKey(PLAINTEXT, K_1);
    addRC(PLAINTEXT, ROUND);

    for (int i=0; i<5; i++){
        ROUND++;
        SBOX(PLAINTEXT);
        MixColumn(PLAINTEXT);
        addRC(PLAINTEXT, ROUND);
        addRoundKey(PLAINTEXT, K_1);
    }
    	SBOX(PLAINTEXT);
    	M_p(PLAINTEXT);
    	INV_SBOX(PLAINTEXT);

    for (int i=0; i<5; i++){
        ROUND++;
        addRoundKey(PLAINTEXT, K_1);
        addRC(PLAINTEXT, ROUND);
        inv_MixColumn(PLAINTEXT);
        INV_SBOX(PLAINTEXT);
    }

    ROUND++;
    addRC(PLAINTEXT, ROUND);
    addRoundKey(PLAINTEXT, K_1);
    KeyExpansion(K_0);
    addRoundKey(PLAINTEXT, K_0);
}












