#include "instruct.h"

Instruct::Instruct(int cyc_deg, int num_primes, int max_prime){
	globalLevel = 0;
	n = new ntru(cyc_deg, num_primes, max_prime, 1, 16);
	c = new myCRT;
}

void Instruct::GenerateParameters(){
		n->Func_SetModulus();
		n->Func_ModulusFindRing(Num_Primes, Max_Prime, Dif_Prime, *n->Pointer_PolyModulus());
		c->SetModulus(*n->Pointer_PolyModulus());
		c->ComputeFactors(n->Pointer_GP()->D_FactDegree, n->Pointer_GP()->FactSize);
		c->CalculateMs();
		c->CalculateNs();
		c->CalculateMxNs();
}

void Instruct::SaveParameters(string in){
	stringstream ss;
	ss << Modulus_M;
	string str = ss.str();
	str = in;

	std::fstream fs;
	fs.open (str.c_str(), std::fstream::out | std::fstream::trunc);

	fs << *n->Pointer_PolyModulus() << endl;

	for(int i=0; i<Num_Primes; i++)
		fs << *n->Pointer_Q(i) << endl;

	fs << c->Returnmodulus() << endl;
	fs << c->Returnsize() << endl;
	vec_ZZ_pX temp;
	temp = c->Returnfactors();

	for(int i=0; i<c->Returnsize(); i++)
		fs << temp[i] << endl;

	temp = c->ReturnM();
	for(int i=0; i<c->Returnsize(); i++)
		fs << temp[i] << endl;

	temp = c->ReturnN();
	for(int i=0; i<c->Returnsize(); i++)
		fs << temp[i] << endl;

	temp = c->ReturnMxN();
	for(int i=0; i<c->Returnsize(); i++)
		fs << temp[i] << endl;

	fs.close();
}

void Instruct::LoadParameters(string in){

	stringstream ss;
	ss << Modulus_M;
	string str = ss.str();
	str = in;

	std::fstream fs;
	fs.open (str.c_str(), std::fstream::in);
	n->IO_ReadModulus(fs);
	n->IO_ReadQs(fs, Num_Primes);
	c->IOReadAll(fs);
	fs.close();
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void Instruct::GenerateKeys(){
	n->Func_ComputeKeysRingRelinFFT(Num_Primes, 1);
}

void Instruct::SaveKeys(string s){
	std::fstream fs;
	fs.open (s.c_str(), std::fstream::out | std::fstream::trunc);
	fs << *n->Pointer_PublicKey() << endl;
	for(int i=0; i<Num_Primes; i++)
		fs << (*n->Pointer_SecretKey(i)) << endl;

	for(int i=0; i<Num_Primes; i++)
		fs << (*n->Pointer_EvalKey(i)) << endl;

	fs.close();
}

void Instruct::LoadKeys(string s){
	std::fstream fs;
	fs.open (s.c_str(), std::fstream::in);
	n->IO_ReadKeys(fs, Num_Primes);
	fs.close();
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


ZZX Instruct::GenRandMess(int size){
	 return n->Func_CreateMessage(size);
}

ntru *Instruct::Pointer_NTRU(){
	return n;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

InstByte Instruct::InttoInstByte(int n){
	InstByte r;

	for(int i=0; i<8; i++)
		r.b[i].bit = (n>>i)%2;

	return r;
}

int	Instruct::InstBytetoInt(InstByte &x){
	int r = 0;
	for(int i=0; i<8; i++)
		r = r + (to_long(coeff(x.b[i].bit,0))<<i);

	return r;
}

InstInt Instruct::InttoInstInt(int n){
	InstInt r;

	for(int i=0; i<32; i++)
		r.b[i].bit = (n>>i)%2;
	return r;
}

int	Instruct::InstInttoInt(InstInt &x){
	int r = 0;
	for(int i=0; i<32; i++)
		r = r + (to_long(coeff(x.b[i].bit,0))<<i);

	return r;
}

InstBit Instruct::C_AND(InstBit &a, InstBit &b, int level){
	InstBit r;
//	r.bit = n->Arith_MulModZZX(a.bit, b.bit, level);
	r.bit = a.bit * b.bit;
	r.lvl = a.lvl;

	return r;
}

InstByte Instruct::C_AND(InstBit &a, InstByte &b, int level){
	InstByte r;
	for(int i=0; i<8; i++)
		r.b[i] = C_AND(a, b.b[i], level);
	r.lvl = a.lvl;
	return r;
}

InstInt Instruct::C_AND(InstBit &a, InstInt &b, int level){
	InstInt r;
	for(int i=0; i<32; i++)
		r.b[i] = C_AND(a, b.b[i], level);
	r.lvl = a.lvl;
	return r;
}

InstBit Instruct::C_XOR(InstBit &a, InstBit &b){
	InstBit r;
	r.bit = a.bit+b.bit;
	r.lvl = a.lvl;

	return r;
}

InstByte Instruct::C_XOR(InstByte &a, InstByte &b){
	InstByte r;

	for(int i=0; i<8; i++)
		r.b[i] = C_XOR(a.b[i], b.b[i]);

	r.lvl = a.lvl;
	return r;
}

InstInt Instruct::C_XOR(InstInt &a, InstInt &b){
	InstInt r;

	for(int i=0; i<32; i++)
		r.b[i] = C_XOR(a.b[i], b.b[i]);

	r.lvl = a.lvl;
	return r;
}

/*
InstBit Instruct::C_XOR(InstBit &a, InstBit &b){

    return (a+b);
}
*/

/*****************************
*** Comparison ( x = y ) -> z
*** ti = (xi + yi)'
*** z = (t0).(t1)...(t7)
*****************************/
InstBit Instruct::isEqual(InstByte &x, InstByte &y){
	InstBit z;
	InstBit t[8];

	int level = EqualizeLevel(x, y);

	// ti = (xi + yi)'
    for(int i=0; i<8; i++)
		t[i] = (x.b[i] + y.b[i]) + 1;

 	/* The ANDs at depth = 1
		(t0).(t4)
		(t1).(t5)
		(t2).(t6)
		(t3).(t7)
	*/
    for(int i=0; i<4; i++){
    	t[i] = C_AND(t[i], t[i+4], level);
    	t[i] = Process1Level(t[i], level);
    }
    level++;


    /* The ANDs at depth = 2
		(t0).(t2)
		(t1).(t3)
	*/
    for(int i=0; i<2; i++){
		t[i] = C_AND(t[i], t[i+2], level);
		t[i] = Process1Level(t[i], level);
    }
    level++;


    /* The AND at depth = 3
		(t0) . (t1)
	*/
    z = C_AND(t[0], t[1], level);
    z = Process1Level(z, level);
    level++;

    return z;
}

/********************************************
*** Less than ( x < y ) -> z
*** ti = (xi + yi)'
*** ki = (xi)'.(yi)
*** si = (ki).(t(i+1)).(t(i+2))...(t7)
*** z = s0 + s1 + ... + s7
*** x[0], y[0] : least significant bit of a, b
*** x[7], y[7] : most significant bit of a, b
*********************************************/

#ifdef OP_isLessThan
InstBit Instruct::isLessThan(InstByte &x, InstByte &y){
	InstBit z;
	InstBit t[8], k[8], d[7], e[5], c[4];

	int level = EqualizeLevel(x, y);


    for(int i=0; i<8; i++)
    	t[i] = (x.b[i] + y.b[i] + 1);		// ti = (xi + yi)'

///////// Depth = 1, #AND = 8
    for(int i=0; i<8; i++){
    	k[i] = x.b[i]+ 1;
    	k[i] = C_AND(k[i], y.b[i], level);
    	k[i] = Process1Level(k[i], level);
    }
    level++;


///////// Depth = 2, #ANDs = 7
    EqualizeLevel(t[7],t[6]);			// t6.t7
    EqualizeLevel(t[5],t[4]);			// t5.t4
    EqualizeLevel(t[3],t[2]);			// t3.t2
    EqualizeLevel(t[1],k[0]);			// t1.k0
    EqualizeLevel(t[3],k[2]);			// t3.k2
    EqualizeLevel(t[5],k[4]);			// t5.k4
    EqualizeLevel(t[7],k[6]);			// t7.k6				=	s6
    d[0] = C_AND(t[7],t[6], level);			// t6.t7
	d[1] = C_AND(t[5],t[4], level);			// t5.t4
	d[2] = C_AND(t[3],t[2], level);			// t3.t2
	d[3] = C_AND(t[1],k[0], level);			// t1.k0
	d[4] = C_AND(t[3],k[2], level);			// t3.k2
	d[5] = C_AND(t[5],k[4], level);			// t5.k4
	d[6] = C_AND(t[7],k[6], level);			// t7.k6				=	s6

	for(int i=0; i<7; i++)
		d[i] = Process1Level(d[i], level);

	level++;


/////////	Depth = 3, #ANDs = 5
	EqualizeLevel(d[0],k[5]);
	EqualizeLevel(d[5],d[0]);
	EqualizeLevel(d[1],d[0]);
	EqualizeLevel(k[1],d[2]);
	EqualizeLevel(d[3],d[2]);

//	e[0] = C_AND(d[0],k[5], level);		// k5.t6.t7				=	s5
//	e[1] = C_AND(d[5],d[0], level);		// k4.t5.t6.t7			=	s4
	e[2] = C_AND(d[1],d[0], level);		// t4.t5.t6.t7
	e[3] = C_AND(k[1],d[2], level);		// k1.t2.t3
	e[4] = C_AND(d[3],d[2], level);		// k0.t1.t2.t3

	for(int i=0; i<5; i++)
		e[i] = Process1Level(e[i], level);

	level++;


	EqualizeLevel(e[4], d[0]);
	EqualizeLevel(e[4], d[5]);
	EqualizeLevel(e[4], k[5]);
	e[0] = C_AND(d[0],k[5], level);
	e[1] = C_AND(d[5],d[0], level);
/////////	Depth = 4, #ANDs = 4
	EqualizeLevel(k[3],e[2]);
	EqualizeLevel(d[4],e[2]);
	EqualizeLevel(e[3],e[2]);
	EqualizeLevel(e[4],e[2]);


	c[0] = C_AND(k[3],e[2], level);		// k3.t4.t5.t6.t7		=	s3
	c[1] = C_AND(d[4],e[2], level);		// k2.t3.t4.t5.t6.t7	=	s2
	c[2] = C_AND(e[3],e[2], level);		// k1.t2...t5.t6.t7		=	s1
	c[3] = C_AND(e[4],e[2], level);		// k0.t1...t5.t6.t7		=	s0

//	for(int i=0; i<4; i++)
//		c[i] = Process1Level(c[i], level);

//	level++;

	// z = s0 + s1 + ... + s7
	EqualizeLevel(c[3], e[0]);
	EqualizeLevel(c[3], e[1]);
	EqualizeLevel(c[3], d[6]);
	EqualizeLevel(c[3], k[7]);

	z = (c[3] + c[2] + c[1] + c[0] + e[1] + e[0] + d[6] + k[7]);

	z  = Process1Level(z, level);
	level++;


    return z;
}
#else
InstBit Instruct::isLessThan(InstByte &x, InstByte &y){
	InstBit z;
	InstBit t[8], k[8], d[7], e[5], c[4];

	int level = EqualizeLevel(x, y);


    for(int i=0; i<8; i++)
    	t[i] = (x.b[i] + y.b[i] + 1);		// ti = (xi + yi)'

///////// Depth = 1, #AND = 8
    for(int i=0; i<8; i++){
    	k[i] = x.b[i]+ 1;
    	k[i] = C_AND(k[i], y.b[i], level);
    	k[i] = Process1Level(k[i], level);
    }
    level++;


///////// Depth = 2, #ANDs = 7
    EqualizeLevel(t[7],t[6]);			// t6.t7
    EqualizeLevel(t[5],t[4]);			// t5.t4
    EqualizeLevel(t[3],t[2]);			// t3.t2
    EqualizeLevel(t[1],k[0]);			// t1.k0
    EqualizeLevel(t[3],k[2]);			// t3.k2
    EqualizeLevel(t[5],k[4]);			// t5.k4
    EqualizeLevel(t[7],k[6]);			// t7.k6				=	s6
    d[0] = C_AND(t[7],t[6], level);			// t6.t7
	d[1] = C_AND(t[5],t[4], level);			// t5.t4
	d[2] = C_AND(t[3],t[2], level);			// t3.t2
	d[3] = C_AND(t[1],k[0], level);			// t1.k0
	d[4] = C_AND(t[3],k[2], level);			// t3.k2
	d[5] = C_AND(t[5],k[4], level);			// t5.k4
	d[6] = C_AND(t[7],k[6], level);			// t7.k6				=	s6

	for(int i=0; i<7; i++)
		d[i] = Process1Level(d[i], level);

	level++;


/////////	Depth = 3, #ANDs = 5
	EqualizeLevel(d[0],k[5]);
	EqualizeLevel(d[5],d[0]);
	EqualizeLevel(d[1],d[0]);
	EqualizeLevel(k[1],d[2]);
	EqualizeLevel(d[3],d[2]);

	e[0] = C_AND(d[0],k[5], level);		// k5.t6.t7				=	s5
	e[1] = C_AND(d[5],d[0], level);		// k4.t5.t6.t7			=	s4
	e[2] = C_AND(d[1],d[0], level);		// t4.t5.t6.t7
	e[3] = C_AND(k[1],d[2], level);		// k1.t2.t3
	e[4] = C_AND(d[3],d[2], level);		// k0.t1.t2.t3

	for(int i=0; i<5; i++)
		e[i] = Process1Level(e[i], level);

	level++;

/////////	Depth = 4, #ANDs = 4
	EqualizeLevel(k[3],e[2]);
	EqualizeLevel(d[4],e[2]);
	EqualizeLevel(e[3],e[2]);
	EqualizeLevel(e[4],e[2]);


	c[0] = C_AND(k[3],e[2], level);		// k3.t4.t5.t6.t7		=	s3
	c[1] = C_AND(d[4],e[2], level);		// k2.t3.t4.t5.t6.t7	=	s2
	c[2] = C_AND(e[3],e[2], level);		// k1.t2...t5.t6.t7		=	s1
	c[3] = C_AND(e[4],e[2], level);		// k0.t1...t5.t6.t7		=	s0

	for(int i=0; i<4; i++)
		c[i] = Process1Level(c[i], level);

	level++;

	// z = s0 + s1 + ... + s7
	EqualizeLevel(c[3], e[0]);
	EqualizeLevel(c[3], e[1]);
	EqualizeLevel(c[3], d[6]);
	EqualizeLevel(c[3], k[7]);

	z = (c[3] + c[2] + c[1] + c[0] + e[1] + e[0] + d[6] + k[7]);

    return z;
}
#endif

int proc=0;
InstBit Instruct::isLessThan(InstInt &x, InstInt &y){
	InstBit z;
	InstBit k[32], t[32];

	int level = EqualizeLevel(x, y);

	for(int i=0; i<32; i++)
		t[i] = (x.b[i] + y.b[i]+1);
// Depth = 1
	//////////////////////////////////////////
	for(int i=0; i<32; i++){
		k[i] = x.b[i]+1;
		k[i] = C_AND(k[i], y.b[i], level);
		if(i == 31)
			Reduce(k[i],level);
		else{
			k[i] = Process1Level(k[i], level);
			proc++;
		}
	}

	for(int i=2; i<32; i+=2){
		t[i] = C_AND(t[i], t[i+1], level);
		t[i] = Process1Level(t[i], level);
		proc++;
	}
	//////////////////////////////////////////
	level++;
	for(int i=0; i<32; i++){
		if(k[i].lvl != 1)
			Update_ValueLevel(k[i], level);
		if(t[i].lvl != 1)
			Update_ValueLevel(t[i], level);
	}
#ifdef __DEBUG__
	////////////DEBUG////////////
									cout << "K:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(k[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
									}
									cout << endl;
									cout << "T:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(t[i], level);
									if(coeff(tt.bit,0) == 0)
										cout << "0";
									else
										cout << "1";
									}
									cout << endl << endl;
		////////////DEBUG////////////
#endif









// Depth = 2
	//////////////////////////////////////////
	for(int i=0; i<32; i+=2){
		k[i] = C_AND(k[i],t[i+1], level);
		if(i == 30)
			Reduce(k[i],level);
		else{
			k[i] = Process1Level(k[i], level);
			proc++;
		}
	}

	for(int i=1; i<32; i+=2){
		EqualizeLevel(k[i-1], k[i]);
		k[i] = (k[i-1]+k[i]);
	}
	for(int i=4; i<32; i+=4){
		t[i] = C_AND(t[i],t[i+2], level);
		t[i] = Process1Level(t[i], level);
		proc++;
	}
	//////////////////////////////////////////
	level++;
	for(int i=0; i<32; i++){
		if(k[i].lvl != 2)
			Update_ValueLevel(k[i], level);
		if(t[i].lvl != 2)
			Update_ValueLevel(t[i], level);
	}
#ifdef __DEBUG__
	////////////DEBUG////////////
									cout << "K:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(k[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
									}
									cout << endl;
									cout << "T:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(t[i], level);
									if(coeff(tt.bit,0) == 0)
										cout << "0";
									else
										cout << "1";
									}
									cout << endl << endl;
		////////////DEBUG////////////
#endif










// Depth = 3, #AND = 8, #Relin = 7
	//////////////////////////////////////////
	for(int i=1; i<32; i+=4){
		k[i] = C_AND(k[i], t[i+1], level);
		if(i == 29)
			Reduce(k[i],level);
		else{
			k[i] = Process1Level(k[i], level);
			proc++;
		}
	}

	for(int i=3; i<32; i+=4){
		EqualizeLevel(k[i], k[i-2]);
		k[i] = (k[i]+k[i-2]);
	}

	for(int i=8; i<32; i+=8){
		t[i] = C_AND(t[i], t[i+4], level);
		t[i] = Process1Level(t[i], level);
		proc++;
	}
	//////////////////////////////////////////
	level++;
	for(int i=0; i<32; i++){
		if(k[i].lvl != 3)
			Update_ValueLevel(k[i], level);
		if(t[i].lvl != 3)
			Update_ValueLevel(t[i], level);
	}

#ifdef __DEBUG__
	////////////DEBUG////////////
									cout << "K:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(k[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
									}
									cout << endl;
									cout << "T:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(t[i], level);
									if(coeff(tt.bit,0) == 0)
										cout << "0";
									else
										cout << "1";
									}
									cout << endl << endl;
		////////////DEBUG////////////
#endif
// Depth = 4, #AND = 4, #Relin = 3
	//////////////////////////////////////////
	for(int i=3; i<32; i+=8){
		k[i] = C_AND(k[i], t[i+1], level);
		if(i == 27)
			Reduce(k[i],level);
		else{
			k[i] = Process1Level(k[i], level);
			proc++;
		}
	}
	for(int i=7; i<32; i+=8){
		EqualizeLevel(k[i], k[i-4]);
		k[i] = (k[i]+k[i-4]);
	}

	t[16] = C_AND(t[16],t[24], level);
	t[16] = Process1Level(t[16], level);
	proc++;
	//////////////////////////////////////////
	level++;
	for(int i=0; i<32; i++){
		if(k[i].lvl != 4)
			Update_ValueLevel(k[i], level);
		if(t[i].lvl != 4)
			Update_ValueLevel(t[i], level);
	}

#ifdef __DEBUG__
	////////////DEBUG////////////
									cout << "K:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(k[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
									}
									cout << endl;
									cout << "T:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(t[i], level);
									if(coeff(tt.bit,0) == 0)
										cout << "0";
									else
										cout << "1";
									}
									cout << endl << endl;
		////////////DEBUG////////////
#endif
// Depth = 5, #AND = 2, #Relin = 1
	//////////////////////////////////////////
		k[7] = C_AND(k[7], t[8], level);
		k[23] = C_AND(k[23], t[24], level);

		EqualizeLevel(k[15], k[7]);
		EqualizeLevel(k[31], k[23]);

		k[15] = (k[15]+k[7]);
		k[31] = (k[31]+k[23]);
	//////////////////////////////////////////

		k[15] = Process1Level(k[15], level);
		proc++;
		level++;
		for(int i=0; i<32; i++){
			if(k[i].lvl != 5)
				Update_ValueLevel(k[i], level);
			if(t[i].lvl != 5)
				Update_ValueLevel(t[i], level);
		}



#ifdef __DEBUG__
		////////////DEBUG////////////
										cout << "K:\t";
										for(int i=0; i<32; i++){
											InstBit tt = Decrypt(k[i], level);
											if(coeff(tt.bit,0) == 0)
												cout << "0";
											else
												cout << "1";
										}
										cout << endl;
										cout << "T:\t";
										for(int i=0; i<32; i++){
											InstBit tt = Decrypt(t[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
										}
										cout << endl << endl;
			////////////DEBUG////////////
#endif




// Depth = 6, #AND = 1
	//////////////////////////////////////////
		k[15] = C_AND(k[15], t[16], level);
	//////////////////////////////////////////

	z = (k[31] + k[15]);
	z = Process1Level(z, level);
	proc++;
	//cout << "my proc:\t"<< proc << endl;
	proc = 0;
#ifdef __DEBUG__
	////////////DEBUG////////////
									cout << "K:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(k[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
									}
									cout << endl;
									cout << "T:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(t[i], level);
									if(coeff(tt.bit,0) == 0)
										cout << "0";
									else
										cout << "1";
									}
									cout << endl << endl;
		////////////DEBUG////////////
#endif
	return z;
}

/*
InstBit Instruct::isLessThan(InstInt &x, InstInt &y){
	InstBit z;
	InstBit k[32], t[32];

	int level = EqualizeLevel(x, y);

	for(int i=0; i<32; i++)
		t[i] = (x.b[i] + y.b[i]+1);
// Depth = 1
	//////////////////////////////////////////
	for(int i=0; i<32; i++){
		k[i] = x.b[i]+1;
		k[i] = C_AND(k[i], y.b[i], level);
		k[i] = Process1Level(k[i], level);
	}

	for(int i=2; i<32; i+=2){
		t[i] = C_AND(t[i], t[i+1], level);
		t[i] = Process1Level(t[i], level);
	}
	//////////////////////////////////////////
	level++;
	for(int i=0; i<32; i++){
		if(k[i].lvl != 1)
			Update_ValueLevel(k[i], level);
		if(t[i].lvl != 1)
			Update_ValueLevel(t[i], level);
	}
#ifdef __DEBUG__
	////////////DEBUG////////////
									cout << "K:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(k[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
									}
									cout << endl;
									cout << "T:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(t[i], level);
									if(coeff(tt.bit,0) == 0)
										cout << "0";
									else
										cout << "1";
									}
									cout << endl << endl;
		////////////DEBUG////////////
#endif









// Depth = 2
	//////////////////////////////////////////
	for(int i=0; i<32; i+=2){
		k[i] = C_AND(k[i],t[i+1], level);
		k[i] = Process1Level(k[i], level);
	}

	for(int i=1; i<32; i+=2){
		EqualizeLevel(k[i-1], k[i]);
		k[i] = (k[i-1]+k[i]);
	}
	for(int i=4; i<32; i+=4){
		t[i] = C_AND(t[i],t[i+2], level);
		t[i] = Process1Level(t[i], level);
	}
	//////////////////////////////////////////
	level++;
	for(int i=0; i<32; i++){
		if(k[i].lvl != 2)
			Update_ValueLevel(k[i], level);
		if(t[i].lvl != 2)
			Update_ValueLevel(t[i], level);
	}
#ifdef __DEBUG__
	////////////DEBUG////////////
									cout << "K:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(k[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
									}
									cout << endl;
									cout << "T:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(t[i], level);
									if(coeff(tt.bit,0) == 0)
										cout << "0";
									else
										cout << "1";
									}
									cout << endl << endl;
		////////////DEBUG////////////
#endif










// Depth = 3, #AND = 8, #Relin = 7
	//////////////////////////////////////////
	for(int i=1; i<32; i+=4){
		k[i] = C_AND(k[i], t[i+1], level);
		k[i] = Process1Level(k[i], level);
	}

	for(int i=3; i<32; i+=4){
		EqualizeLevel(k[i], k[i-2]);
		k[i] = (k[i]+k[i-2]);
	}

	for(int i=8; i<32; i+=8){
		t[i] = C_AND(t[i], t[i+4], level);
		t[i] = Process1Level(t[i], level);
	}
	//////////////////////////////////////////
	level++;
	for(int i=0; i<32; i++){
		if(k[i].lvl != 3)
			Update_ValueLevel(k[i], level);
		if(t[i].lvl != 3)
			Update_ValueLevel(t[i], level);
	}

#ifdef __DEBUG__
	////////////DEBUG////////////
									cout << "K:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(k[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
									}
									cout << endl;
									cout << "T:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(t[i], level);
									if(coeff(tt.bit,0) == 0)
										cout << "0";
									else
										cout << "1";
									}
									cout << endl << endl;
		////////////DEBUG////////////
#endif
// Depth = 4, #AND = 4, #Relin = 3
	//////////////////////////////////////////
	for(int i=3; i<32; i+=8){
		k[i] = C_AND(k[i], t[i+1], level);
		k[i] = Process1Level(k[i], level);
	}
	for(int i=7; i<32; i+=8){
		EqualizeLevel(k[i], k[i-4]);
		k[i] = (k[i]+k[i-4]);
	}

	t[16] = C_AND(t[16],t[24], level);
	t[16] = Process1Level(t[16], level);
	//////////////////////////////////////////
	level++;
	for(int i=0; i<32; i++){
		if(k[i].lvl != 4)
			Update_ValueLevel(k[i], level);
		if(t[i].lvl != 4)
			Update_ValueLevel(t[i], level);
	}

#ifdef __DEBUG__
	////////////DEBUG////////////
									cout << "K:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(k[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
									}
									cout << endl;
									cout << "T:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(t[i], level);
									if(coeff(tt.bit,0) == 0)
										cout << "0";
									else
										cout << "1";
									}
									cout << endl << endl;
		////////////DEBUG////////////
#endif
// Depth = 5, #AND = 2, #Relin = 1
	//////////////////////////////////////////
		k[7] = C_AND(k[7], t[8], level);
		k[23] = C_AND(k[23], t[24], level);

		EqualizeLevel(k[15], k[7]);
		EqualizeLevel(k[31], k[23]);

		k[15] = (k[15]+k[7]);
		k[31] = (k[31]+k[23]);
	//////////////////////////////////////////

		k[15] = Process1Level(k[15], level);

		level++;
		for(int i=0; i<32; i++){
			if(k[i].lvl != 5)
				Update_ValueLevel(k[i], level);
			if(t[i].lvl != 5)
				Update_ValueLevel(t[i], level);
		}



#ifdef __DEBUG__
		////////////DEBUG////////////
										cout << "K:\t";
										for(int i=0; i<32; i++){
											InstBit tt = Decrypt(k[i], level);
											if(coeff(tt.bit,0) == 0)
												cout << "0";
											else
												cout << "1";
										}
										cout << endl;
										cout << "T:\t";
										for(int i=0; i<32; i++){
											InstBit tt = Decrypt(t[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
										}
										cout << endl << endl;
			////////////DEBUG////////////
#endif




// Depth = 6, #AND = 1
	//////////////////////////////////////////
		k[15] = C_AND(k[15], t[16], level);
	//////////////////////////////////////////

	z = (k[31] + k[15]);
	z = Process1Level(z, level);

#ifdef __DEBUG__
	////////////DEBUG////////////
									cout << "K:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(k[i], level);
										if(coeff(tt.bit,0) == 0)
											cout << "0";
										else
											cout << "1";
									}
									cout << endl;
									cout << "T:\t";
									for(int i=0; i<32; i++){
										InstBit tt = Decrypt(t[i], level);
									if(coeff(tt.bit,0) == 0)
										cout << "0";
									else
										cout << "1";
									}
									cout << endl << endl;
		////////////DEBUG////////////
#endif
	return z;
}
*/





InstByte Instruct::Add(InstByte &x, InstByte &y){
	InstByte r;
	InstBit p[8], p1[6], p2[4], g[8], pg[7], pg1[6], pg2[5], pg3[9], pg4, c[9];
	InstBit Sum;

	int level = EqualizeLevel(x, y);

	for(int i=0; i<8; i++)
		p[i] = (x.b[i] + y.b[i]);			// pi = (xi + yi)

//////// Depth = 1, #AND = 8
	for(int i=0; i<8; i++){
		g[i] = C_AND(x.b[i],y.b[i], level);		// gi = (xi).(yi)
		g[i] = Process1Level(g[i], level);
	}

	// Depth = 1, #AND = 6
	for(int i=0; i<6; i++){
		p1[i] = C_AND(p[i+1], p[i+2], level);
		p1[i] = Process1Level(p1[i], level);
	}
	level++;
										// p1[0] = P2.P1
										// p1[1] = P3.P2
										// p1[2] = P4.P3
										// p1[3] = P5.P4
										// p1[4] = P6.P5
										// p1[5] = P7.P6




//////// Depth = 2, #AND = 7
	for(int i=0; i<7; i++){
		EqualizeLevel(p[i+1],g[i]);
		pg[i] = C_AND(p[i+1], g[i], level);
		pg[i] = Process1Level(pg[i], level);
	}
										// pg[0] = P1.G0
										// pg[1] = P2.G1
										// pg[2] = P3.G2
										// pg[3] = P4.G3
										// pg[4] = P5.G4
										// pg[5] = P6.G5
										// pg[6] = P7.G6

	// Depth = 2, #AND = 4
	for(int i=0; i<4; i++){
		EqualizeLevel(p1[i],p1[i+2]);
		p2[i] = C_AND(p1[i],p1[i+2], level);
		p2[i] = Process1Level(p2[i], level);
	}
										// p2[0] = P4.P3.P2.P1
										// p2[1] = P5.P4.P3.P2
										// p2[2] = P6.P5.P4.P3
										// p2[3] = P7.P6.P5.P4
	// Depth = 2, #AND = 6
	for(int i=0; i<6; i++){
		EqualizeLevel(p1[i], g[i]);
		pg1[i] = C_AND(p1[i], g[i], level);
		pg1[i] = Process1Level(pg1[i], level);
	}
	level++;
										// pg1[0] = P2.P1.G0
										// pg1[1] = P3.P2.G1
										// pg1[2] = P4.P3.G2
										// pg1[3] = P5.P4.G3
										// pg1[4] = P6.P5.G4
										// pg1[5] = P7.P6.G5


//////// Depth = 3, #AND = 5
	for(int i=0; i<5; i++){
		EqualizeLevel(p1[i+1], pg[i]);
		pg2[i] = C_AND(p1[i+1], pg[i], level);
		pg2[i] = Process1Level(pg2[i], level);
	}
										// pg2[0] = P3.P2.P1.G0
										// pg2[1] = P4.P3.P2.G1
										// pg2[2] = P5.P4.P3.G2
										// pg2[3] = P6.P5.P4.G3
										// pg2[4] = P7.P6.P5.G4
	// Depth = 3, #AND = 4
	for(int i=0; i<4; i++){
		EqualizeLevel(p2[i],g[i]);
		pg3[i] = C_AND(p2[i],g[i], level);
		pg3[i] = Process1Level(pg3[i], level);
	}
										// pg3[0] = P4.P3.P2.P1.G0
										// pg3[1] = P5.P4.P3.P2.G1
										// pg3[2] = P6.P5.P4.P3.G2
										// pg3[3] = P7.P6.P5.P4.G3
	// Depth = 3, #AND = 3
	for(int i=4; i<7; i++){
		EqualizeLevel(p2[i-3],pg[i-4]);
		pg3[i] = C_AND(p2[i-3],pg[i-4], level);
		pg3[i] = Process1Level(pg3[i], level);
	}
											// pg3[4] = P5.P4.P3.P2.P1.G0
											// pg3[5] = P6.P5.P4.P3.P2.G1
											// pg3[6] = P7.P6.P5.P4.P3.G2
	// Depth = 3, #AND = 2
	for(int i=7; i<9; i++){
		EqualizeLevel(p2[i-5], pg1[i-7]);
		pg3[i] = C_AND(p2[i-5], pg1[i-7], level);
		pg3[i] = Process1Level(pg3[i], level);
	}
	level++;
											// pg3[7] = P6.P5.P4.P3.P2.P1.G0
											// pg3[8] = P7.P6.P5.P4.P3.P2.G1
//////// Depth = 4, #AND = 1
	//pg4 = C_AND(p2[3],pg2[0]);				// pg4 = P7.P6.P5.P4.P3.P2.P1.G0
	// Finally calculate Ci
	c[0].bit = 0;	// First Carry-in is 0 for addition
	c[0].lvl = 0;
	EqualizeLevel(pg3[0], c[0]);

	EqualizeLevel(pg3[0], g[0]);
	c[1] = g[0];

	EqualizeLevel(pg3[0], g[1]);
	EqualizeLevel(pg3[0], pg[0]);
	c[2] = (g[1] + pg[0]);

	EqualizeLevel(pg3[0], g[2]);
	EqualizeLevel(pg3[0], pg[1]);
	EqualizeLevel(pg3[0], pg1[0]);
	c[3] = (g[2] + pg[1] + pg1[0]);

	EqualizeLevel(pg3[0], g[3]);
	EqualizeLevel(pg3[0], pg[2]);
	EqualizeLevel(pg3[0], pg1[1]);
	EqualizeLevel(pg3[0], pg2[0]);
	c[4] = (g[3] + pg[2] + pg1[1] + pg2[0]);

	EqualizeLevel(pg3[0], g[4]);
	EqualizeLevel(pg3[0], pg[3]);
	EqualizeLevel(pg3[0], pg1[2]);
	EqualizeLevel(pg3[0], pg2[1]);
	EqualizeLevel(pg3[0], pg3[0]);
	c[5] = (g[4] + pg[3] + pg1[2] + pg2[1] + pg3[0]);

	EqualizeLevel(pg3[0], g[5]);
	EqualizeLevel(pg3[0], pg[4]);
	EqualizeLevel(pg3[0], pg1[3]);
	EqualizeLevel(pg3[0], pg2[2]);
	EqualizeLevel(pg3[0], pg3[1]);
	EqualizeLevel(pg3[0], pg3[4]);
	c[6] = (g[5] + pg[4] + pg1[3] + pg2[2] + pg3[1] + pg3[4]);

	EqualizeLevel(pg3[0], g[6]);
	EqualizeLevel(pg3[0], pg[5]);
	EqualizeLevel(pg3[0], pg1[4]);
	EqualizeLevel(pg3[0], pg2[3]);
	EqualizeLevel(pg3[0], pg3[2]);
	EqualizeLevel(pg3[0], pg3[5]);
	EqualizeLevel(pg3[0], pg3[7]);
	c[7] = (g[6] + pg[5] + pg1[4] + pg2[3] + pg3[2] + pg3[5] + pg3[7]);
	//c[8] = (g[7] + pg[6] + pg1[5] + pg2[4] + pg3[3] + pg3[6] + pg3[8] + pg4)%2;	// The carry-out
	// Si = Pi + Ci

	for(int i=0; i<8; i++){
		EqualizeLevel(p[i], c[i]);
		r.b[i] = (p[i] + c[i]);
	}
	return r;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


void Instruct::GreedySort(InstByte *out, InstByte *in, int size){

	if(size == 0)
		return;
	else if(size == 1)
		out[0] = in[0];
	else if(size == 2){


		InstBit A = isLessThan(in[0], in[1]);

        Update_GlobalLevel(level_isLessThan_Byte);
        Update_TableLevel();
        Update_ValueLevel(in[0]);
        Update_ValueLevel(in[1]);

		InstBit B = A + 1;

		//////////////////////////////////

		out[0] = C_AND(A, in[0], globalLevel) + C_AND(B, in[1], globalLevel);//t[0] + t[1];
		out[1] = C_AND(A, in[1], globalLevel) + C_AND(B, in[0], globalLevel); //t[2] + t[3];

		out[0] = Process1Level(out[0], globalLevel);
		out[1] = Process1Level(out[1], globalLevel);

		Update_GlobalLevel(level_and);
		Update_TableLevel();
		Update_ValueLevel(out[0]);
		Update_ValueLevel(out[1]);
	}
	else{

		InstBit** M = new InstBit*[size];
		for(int i=0; i<size; i++)
			M[i] = new InstBit[size];

		for(int i=0; i<size; i++){
			M[i][i].bit = 0;
			M[i][i].lvl = level_isLessThan_Byte;
			for(int j=i+1; j<size; j++){
				M[i][j] = isLessThan(in[i], in[j]);		// ( x[i] < x[j] )
				M[j][i] = M[i][j] + 1;				// ( x[i] < x[j] )'
			}
		}
		Update_GlobalLevel(level_isLessThan_Byte);
		Update_TableLevel();

		int iter = 0;
		int C_row_size = 1;
		int row_cnt = 3;
		for(int i=1; i<size-1; i*=2)
			iter++;
		for(int i=0; i<iter; i++)
			C_row_size *=2;
		C_row_size+=1;

		InstBit ***C = new InstBit**[size];
		for(int i=0; i<size; i++){
			C[i] = new InstBit*[C_row_size];
			for(int j=0; j<C_row_size; j++)
				C[i][j] = new InstBit[size/2];
		}


		////// Depth 4
		for(int i=0; i<size; i++){
			int j1 = (i+1)%size;
			int end = (i-1+size)%size;

			while(j1!=end && j1!=i){

				int j2 = (j1+1)%size;
				int k = ((j1-i+size)%size)/2;

				InstBit mytemp = C_AND(M[i][j1], M[i][j2], globalLevel);
				InstBit mytemp2 = (M[i][j1]+M[i][j2]);

				mytemp = Process1Level(mytemp, globalLevel);
				EqualizeLevel(mytemp, mytemp2);

				C[i][0][k] = mytemp;
				C[i][1][k] = mytemp2;
				C[i][2][k] = (mytemp + mytemp2 + 1);

				j1 = (j2+1)%size;
				// RELIN HERE

			}
			if(j1==end){
				int k = size/2-1;
				C[i][0][k] 		= M[i][j1];
				C[i][1][k] 		= M[j1][i];
				C[i][2][k].bit 	= 0;
				C[i][2][k].lvl 	= 0;
			}
		}
		Update_GlobalLevel(level_and);
		Update_TableLevel();

		int new_row_cnt = (row_cnt-1)*2 + 1;
		int new_size = size/2;

		////// Depth 5
		for(int k=new_size; k>1; k=(k+1)/2){		// exactly lg(size-1/2) times do this
			for(int i=0; i<size; i++){			// do this for every x
				for(int l=0; l<k; l+=2){			// do this for every consecutive pair

					if(l == k-1){
						for(int a=0; a<row_cnt; a++)
							C[i][a][l/2] = C[i][a][l];
						for(int a=row_cnt; a<new_row_cnt; a++){
							C[i][a][l/2].bit = 0;
							C[i][a][l/2].lvl = globalLevel;
						}
					}
					else{
						InstBit *temp = new InstBit[new_row_cnt];
						for(int a=0; a<new_row_cnt; a++){
							temp[a].bit = 0;
							temp[a].lvl = 0;
						}

						for(int a=0; a<row_cnt; a++){	// do this for every row of C
							for(int b=0; b<row_cnt; b++){
								temp[a+b] = temp[a+b] + C_AND(C[i][a][l], C[i][b][l+1], globalLevel);
							// RELIN HERE
							}

						}

						// Set C according to temp
						for(int a=0; a<new_row_cnt; a++){
							temp[a] = Process1Level(temp[a], globalLevel);
							C[i][a][l/2] = temp[a];
						}

						delete[] temp;
					}
				}
			}
			row_cnt = new_row_cnt;
			new_row_cnt = (new_row_cnt-1)*2 + 1;
			Update_GlobalLevel(level_and);
			Update_TableLevel();
		}

		////// Depth 5 + lg(size-1/2)
		////// Depth 4 + lg(size-1)
		InstByte **Btemp = new InstByte*[size];
		for(int i=0; i<size; i++)
			Btemp[i] = new InstByte[size];

		for(int i=0; i<size; i++){
			for(int j=0; j<size; j++){
				Btemp[i][j] = C_AND(C[i][j][0], in[i], globalLevel);
			}
		}

		////// Depth 5 + lg(size-1)
	for(int i=0; i<size; i++){
			out[i].b[0].bit = 0;
			out[i].b[1].bit = 0;
			out[i].b[2].bit = 0;
			out[i].b[3].bit = 0;

			out[i].b[4].bit = 0;
			out[i].b[5].bit = 0;
			out[i].b[6].bit = 0;
			out[i].b[7].bit = 0;

			for(int j=0; j<size; j++){
				out[i] = out[i] + Btemp[j][i];
						// RELIN HERE
			}
			out[i] = Process1Level(out[i], globalLevel);
		}
		Update_GlobalLevel(level_and);
		Update_TableLevel();


		for(int i=0; i<size; i++)
			delete []Btemp[i];
		delete [] Btemp;

		for(int i=0; i<size;i++){
			delete[] M[i];
			for(int j=0; j<C_row_size; j++)
				delete[] C[i][j];
			delete[] C[i];
		}
		delete[] M;
		delete[] C;
	}
}

void Instruct::GreedySort(InstInt *out, InstInt *in, int size){

	if(size == 0)
		return;
	else if(size == 1)
		out[0] = in[0];
	else if(size == 2){

		InstBit A = isLessThan(in[0], in[1]);

        Update_GlobalLevel(level_isLessThan_Int);
        Update_TableLevel();
        Update_ValueLevel(in[0]);
        Update_ValueLevel(in[1]);

		InstBit B = A + 1;

		//////////////////////////////////
		out[0] = C_AND(A, in[0], globalLevel) + C_AND(B, in[1], globalLevel);//t[0] + t[1];
		out[1] = C_AND(A, in[1], globalLevel) + C_AND(B, in[0], globalLevel); //t[2] + t[3];

		out[0] = Process1Level(out[0], globalLevel);
		out[1] = Process1Level(out[1], globalLevel);

		Update_GlobalLevel(level_and);
		Update_TableLevel();
		Update_ValueLevel(out[0]);
		Update_ValueLevel(out[1]);
	}
	else{
		InstBit** M = new InstBit*[size];
		for(int i=0; i<size; i++)
			M[i] = new InstBit[size];

		for(int i=0; i<size; i++){
			M[i][i].bit = 0;
			M[i][i].lvl = level_isLessThan_Int;
			for(int j=i+1; j<size; j++){
				M[i][j] = isLessThan(in[i], in[j]);		// ( x[i] < x[j] )
				M[j][i] = M[i][j] + 1;				// ( x[i] < x[j] )'
			}
		}
		Update_GlobalLevel(level_isLessThan_Int);
		Update_TableLevel();

		int iter = 0;
		int C_row_size = 1;
		int row_cnt = 3;
		for(int i=1; i<size-1; i*=2)
			iter++;
		for(int i=0; i<iter; i++)
			C_row_size *=2;
		C_row_size+=1;


		InstBit ***C = new InstBit**[size];
		for(int i=0; i<size; i++){
			C[i] = new InstBit*[C_row_size];
			for(int j=0; j<C_row_size; j++)
				C[i][j] = new InstBit[size/2];
		}

		////// Depth 4
		for(int i=0; i<size; i++){
			int j1 = (i+1)%size;
			int end = (i-1+size)%size;

			while(j1!=end && j1!=i){

				int j2 = (j1+1)%size;
				int k = ((j1-i+size)%size)/2;

				InstBit mytemp = C_AND(M[i][j1], M[i][j2], globalLevel);
				InstBit mytemp2 = (M[i][j1]+M[i][j2]);

				mytemp = Process1Level(mytemp, globalLevel);
				EqualizeLevel(mytemp, mytemp2);

				C[i][0][k] = mytemp;
				C[i][1][k] = mytemp2;
				C[i][2][k] = (mytemp + mytemp2 + 1);

				j1 = (j2+1)%size;
				// RELIN HERE

			}
			if(j1==end){
				int k = size/2-1;
				C[i][0][k] 		= M[i][j1];
				C[i][1][k] 		= M[j1][i];
				C[i][2][k].bit 	= 0;
				C[i][2][k].lvl 	= 0;
			}
		}
		Update_GlobalLevel(level_and);
		Update_TableLevel();

		int new_row_cnt = (row_cnt-1)*2 + 1;
		int new_size = size/2;

		////// Depth 5
		for(int k=new_size; k>1; k=(k+1)/2){		// exactly lg(size-1/2) times do this
			for(int i=0; i<size; i++){			// do this for every x
				for(int l=0; l<k; l+=2){			// do this for every consecutive pair

					if(l == k-1){
						for(int a=0; a<row_cnt; a++)
							C[i][a][l/2] = C[i][a][l];
						for(int a=row_cnt; a<new_row_cnt; a++){
							C[i][a][l/2].bit = 0;
							C[i][a][l/2].lvl = globalLevel;
						}
					}
					else{
						InstBit *temp = new InstBit[new_row_cnt];
						for(int a=0; a<new_row_cnt; a++){
							temp[a].bit = 0;
							temp[a].lvl = 0;
						}

						for(int a=0; a<row_cnt; a++){	// do this for every row of C
							for(int b=0; b<row_cnt; b++){
								temp[a+b] = temp[a+b] + C_AND(C[i][a][l], C[i][b][l+1], globalLevel);
							// RELIN HERE
							}

						}

						// Set C according to temp
						for(int a=0; a<new_row_cnt; a++){
							temp[a] = Process1Level(temp[a], globalLevel);
							C[i][a][l/2] = temp[a];
						}

						delete[] temp;
					}
				}
			}
			row_cnt = new_row_cnt;
			new_row_cnt = (new_row_cnt-1)*2 + 1;
			Update_GlobalLevel(level_and);
			Update_TableLevel();
		}

		////// Depth 5 + lg(size-1/2)
		////// Depth 4 + lg(size-1)
		InstInt **Btemp = new InstInt*[size];
		for(int i=0; i<size; i++)
			Btemp[i] = new InstInt[size];

		for(int i=0; i<size; i++){
			for(int j=0; j<size; j++){
				Btemp[i][j] = C_AND(C[i][j][0], in[i], globalLevel);
			}
		}


		////// Depth 5 + lg(size-1)
        for(int i=0; i<size; i++){
			out[i].b[0].bit = 0;
			out[i].b[1].bit = 0;
			out[i].b[2].bit = 0;
			out[i].b[3].bit = 0;

			out[i].b[4].bit = 0;
			out[i].b[5].bit = 0;
			out[i].b[6].bit = 0;
			out[i].b[7].bit = 0;

			for(int j=0; j<size; j++){
				out[i] = out[i] + Btemp[j][i];
						// RELIN HERE
			}
			out[i] = Process1Level(out[i], globalLevel);
		}
		Update_GlobalLevel(level_and);
		Update_TableLevel();

		for(int i=0; i<size; i++)
			delete []Btemp[i];
		delete [] Btemp;

		for(int i=0; i<size;i++){
			delete[] M[i];
			for(int j=0; j<C_row_size; j++)
				delete[] C[i][j];
			delete[] C[i];
		}
		delete[] M;
		delete[] C;
	}

}


/////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

InstByte Instruct::HammingWeight(InstBit *A, int size){

	InstByte result;
	int level = globalLevel;

	if(size==4)
	{
		InstBit s1 = A[0] + A[1] + A[2];

		result.b[0] = s1 + A[3];

        // Depth 1
		InstBit c1 = C_AND(A[0], A[1], level) + C_AND(A[0], A[2], level) + C_AND(A[1], A[2], level);
        InstBit c2 = C_AND(s1, A[3], level);

		result.b[1] = c1 + c2;

		// Finalize
		result.b[1] = Process1Level(result.b[1], level);
		level++;
		Update_ValueLevel(result.b[0], level);

		result.lvl = level;
	}
	if(size==8)
	{
		InstBit s1 = A[0] + A[1] + A[2];
		InstBit s2 = A[3] + A[4] + A[5];
		InstBit s3 = A[6] + A[7];
        InstBit s11 = s1 + s2 + s3;

		result.b[0] = s11;

        // Depth 1
		InstBit c1 = C_AND(A[0], A[1], level) + C_AND(A[0], A[2], level) + C_AND(A[1], A[2], level);
		c1 = Process1Level(c1, level);
		InstBit c2 = C_AND(A[3], A[4], level) + C_AND(A[3], A[5], level) + C_AND(A[4], A[5], level);
		c2 = Process1Level(c2, level);
		InstBit c3 = C_AND(A[6], A[7], level);
        c3 = Process1Level(c3, level);
        InstBit s21 = c1 + c2 + c3;

        InstBit c11 = C_AND(s1, s2, level) + C_AND(s1, s3, level) + C_AND(s2, s3, level);
        c11 = Process1Level(c11, level);
        InstBit s22 = s21 + c11;

        result.b[1] = s22;
        level++;

		// Depth 2
		InstBit c21 = C_AND(c1, c2, level) + C_AND(c1, c3, level) + C_AND(c2, c3, level);
		InstBit c22 = C_AND(s21, c11, level);
		InstBit s33 = c21 + c22;

		result.b[2] = s33;
        result.b[2] = Process1Level(result.b[2], level);
        level++;
        Update_ValueLevel(result.b[1], level);
		Update_ValueLevel(result.b[0], level);

		result.lvl = level;
	}

    else if(size==16)
	{
		InstBit S[6];
		InstBit C[6];

        // Depth 1
		for(int i=0; i<4; i++)
		{
			S[i] = (A[3*i] + A[3*i+1] + A[3*i+2]);
			C[i] = C_AND(A[3*i], A[3*i+1], level)
                 + C_AND(A[3*i], A[3*i+2], level) + C_AND(A[3*i+1], A[3*i+2], level);
            C[i] = Process1Level(C[i], level);
		}
		S[4] = (A[12] + A[13]);
		C[4] = C_AND(A[12], A[13], level);
		C[4] = Process1Level(C[4], level);

        S[5] = (A[14] + A[15]);
		C[5] = C_AND(A[14], A[15], level);
		C[5] = Process1Level(C[5], level);

        InstBit S1 = (S[0] + S[1] + S[2]);
		InstBit S2 = (S[3] + S[4] + S[5]);
		InstBit S3 = (C[0] + C[1] + C[2]);
		InstBit S4 = (C[3] + C[4] + C[5]);


		InstBit C1 = C_AND(S[0], S[1], level)
                    + C_AND(S[0], S[2], level) + C_AND(S[1], S[2], level);
		C1 = Process1Level(C1, level);

		S[0] = (S3 + S4 + C1);

		S[1] = C_AND(S[3], S[4], level)
                    + C_AND(S[3], S[5], level) + C_AND(S[4], S[5], level);
		S[1] = Process1Level(S[1], level);

		S[2] = C_AND(S1, S2, level);
		S[2] = Process1Level(S[2], level);

		result.b[0] = S1 + S2;

		result.b[1] = S[0] + S[1] + S[2];

		level++;

		// Depth 2
        C[0] = C_AND(C[0], C[1], level)
                    + C_AND(C[0], C[2], level) + C_AND(C[1], C[2], level);
		C[0] = Process1Level(C[0], level);

		C[1] = C_AND(C[3], C[4], level)
                    + C_AND(C[3], C[5], level) + C_AND(C[4], C[5], level);
		C[1] = Process1Level(C[1], level);

		C[2] = C_AND(S3, S4, level)
                    + C_AND(S3, C1, level) + C_AND(S4, C1, level);
		C[2] = Process1Level(C[2], level);

		C[3] = C_AND(S[0], S[1], level)
                    + C_AND(S[0], S[2], level) + C_AND(S[1], S[2], level);
		C[3] = Process1Level(C[3], level);

		S1 = C[0] + C[1];
		S2 = C[2] + C[3];

		result.b[2] = S1 + S2;

		level++;

		// Depth 3

		C[0] = C_AND(C[0], C[1], level);
		C[1] = C_AND(C[2], C[3], level);
		C[2] = C_AND(S1, S2, level);

		result.b[3] = C[0] + C[1] + C[2];
		result.b[3] = Process1Level(result.b[3], level);

		level++;

		Update_ValueLevel(result.b[2], level);
        Update_ValueLevel(result.b[1], level);
		Update_ValueLevel(result.b[0], level);

        result.lvl = level;
	}
	else if(size==32)
	{
		InstBit S[12], C[10];

		// Depth 1
		for(int i=0; i<30; i+=3)
		{
            int j = i/3;
			S[j] = (A[i] + A[i+1] + A[i+2]);
			C[j] = C_AND(A[i], A[i+1], level)
                 + C_AND(A[i], A[i+2], level) + C_AND(A[i+1], A[i+2], level);
            C[j] = Process1Level(C[j], level);
		}
		S[10] = A[30];
		S[11] = A[31];

		InstBit S2[4];
		InstBit C2[8];
		C2[0] = C[9];

		for(int i=0; i<12; i+=3)
		{
            int j = i/3;
			S2[j] = (S[i] + S[i+1] + S[i+2]);
			C2[j+1] = C_AND(S[i], S[i+1], level)
                  + C_AND(S[i], S[i+2], level) + C_AND(S[i+1], S[i+2], level);
            C2[j+1] = Process1Level(C2[j+1], level);
		}

		S[0] = S2[0] + S2[1];
		S[1] = S2[2] + S2[3];

		result.b[0] = S[0] + S[1];

		C2[5] = C_AND(S2[0], S2[1], level);
		C2[5] = Process1Level(C2[5], level);
		C2[6] = C_AND(S2[2], S2[3], level);
		C2[6] = Process1Level(C2[6], level);
		C2[7] = C_AND(S[0], S[1], level);
		C2[7] = Process1Level(C2[7], level);

		level++;

		// Depth 2
		for(int i=0; i<9; i+=3)
		{
            int j = i/3;
			S[j] = (C[i] + C[i+1] + C[i+2]);
			C[j] = C_AND(C[i], C[i+1], level)
                 + C_AND(C[i], C[i+2], level) + C_AND(C[i+1], C[i+2], level);
            C[j] = Process1Level(C[j], level);
		}
		for(int i=0; i<6; i+=3)
		{
            int j = 3 + i/3;
			S[j] = (C2[i] + C2[i+1] + C2[i+2]);
			C[j] = C_AND(C2[i], C2[i+1], level)
                 + C_AND(C2[i], C2[i+2], level) + C_AND(C2[i+1], C2[i+2], level);
            C[j] = Process1Level(C[j], level);
		}
        S[5] = C2[6] + C2[7];
        C[5] = C_AND(C2[6], C2[7], level);
        C[5] = Process1Level(C[5], level);

        S2[0] = (S[0] + S[1] + S[2]);
        S2[1] = (S[3] + S[4] + S[5]);

        result.b[1] = S2[0] + S2[1];

        S[0] = C_AND(S[0], S[1], level)
             + C_AND(S[0], S[2], level) + C_AND(S[1], S[2], level);
        S[0] = Process1Level(S[0], level);
        S[1] = C_AND(S[3], S[4], level)
             + C_AND(S[3], S[5], level) + C_AND(S[4], S[5], level);
        S[1] = Process1Level(S[1], level);

        S[2] = (C[0] + C[1] + C[2]);
        S[3] = (C[3] + C[4] + C[5]);

        S[4] = C_AND(S2[0], S2[1], level);
        S[4] = Process1Level(S[4], level);

        level++;

        // Depth 3
        C[0] = C_AND(C[0], C[1], level)
             + C_AND(C[0], C[2], level) + C_AND(C[1], C[2], level);
        C[0] = Process1Level(C[0], level);
        C[1] = C_AND(C[3], C[4], level)
             + C_AND(C[3], C[5], level) + C_AND(C[4], C[5], level);
        C[1] = Process1Level(C[1], level);
        C[2] = C_AND(S[0], S[1], level)
             + C_AND(S[0], S[2], level) + C_AND(S[1], S[2], level);
        C[2] = Process1Level(C[2], level);
        C[3] = C_AND(S[3], S[4], level);
        C[3] = Process1Level(C[3], level);

        S[0] = (S[0] + S[1] + S[2]);
        S[1] = (S[3] + S[4]);

        result.b[2] = S[0] + S[1];

        C[4] = C_AND(S[0], S[1], level);
        C[4] = Process1Level(C[4], level);

        S[0] = (C[0] + C[1] + C[2]);
        S[1] = (C[3] + C[4]);

        result.b[3] = S[0] + S[1];

        level++;

        // Depth 4
        C[0] = C_AND(C[0], C[1], level)
             + C_AND(C[0], C[2], level) + C_AND(C[1], C[2], level);
        C[1] = C_AND(C[3], C[4], level);
        C[2] = C_AND(S[0], S[1], level);

        result.b[4] = C[0] + C[1] + C[2];
        result.b[4] = Process1Level(result.b[4], level);

        level++;

        Update_ValueLevel(result.b[3], level);
        Update_ValueLevel(result.b[2], level);
        Update_ValueLevel(result.b[1], level);
		Update_ValueLevel(result.b[0], level);

        result.lvl = level;
    }

    else if(size==64)
	{
		InstBit S[22];
		InstBit C[21];

		// Depth 1
		for(int i=0; i<63; i+=3)
		{
            int j = i/3;
			S[j] = (A[i] + A[i+1] + A[i+2]);
			C[j] = C_AND(A[i], A[i+1], level)
                 + C_AND(A[i], A[i+2], level) + C_AND(A[i+1], A[i+2], level);
            C[j] = Process1Level(C[j], level);
		}
		S[21] = A[63];

		InstBit C2[11];

		for(int i=0; i<21; i+=3)
		{
            int j = i/3;
			C2[j] = C_AND(S[i], S[i+1], level)
                  + C_AND(S[i], S[i+2], level) + C_AND(S[i+1], S[i+2], level);
            C2[j] = Process1Level(C2[j], level);
            S[j] = (S[i] + S[i+1] + S[i+2]);
		}
		S[7] = S[21];

        C2[7] = C_AND(S[0], S[1], level)
              + C_AND(S[0], S[2], level) + C_AND(S[1], S[2], level);
        C2[7] = Process1Level(C2[7], level);
        C2[8] = C_AND(S[3], S[4], level)
              + C_AND(S[3], S[5], level) + C_AND(S[4], S[5], level);
        C2[8] = Process1Level(C2[8], level);
        C2[9] = C_AND(S[6], S[7], level);
        C2[9] = Process1Level(C2[9], level);

        S[0] = (S[0] + S[1] + S[2]);
        S[1] = (S[3] + S[4] + S[5]);
        S[2] = (S[6] + S[7]);

        C2[10] = C_AND(S[0], S[1], level)
              + C_AND(S[0], S[2], level) + C_AND(S[1], S[2], level);
        C2[10] = Process1Level(C2[10], level);

        result.b[0] = S[0] + S[1] + S[2];

		level++;

        // Depth 2
		for(int i=0; i<21; i+=3)
		{
            int j = i/3;
            S[j] = (C[i] + C[i+1] + C[i+2]);
			C[j] = C_AND(C[i], C[i+1], level)
                  + C_AND(C[i], C[i+2], level) + C_AND(C[i+1], C[i+2], level);
            C[j] = Process1Level(C[j], level);

		}
		S[7] = (C2[0] + C2[1] + C2[2]);
		C[7] = C_AND(C2[0], C2[1], level)
              + C_AND(C2[0], C2[2], level) + C_AND(C2[1], C2[2], level);
        C[7] = Process1Level(C[7], level);

        S[8] = (C2[3] + C2[4] + C2[5]);
		C[8] = C_AND(C2[3], C2[4], level)
              + C_AND(C2[3], C2[5], level) + C_AND(C2[4], C2[5], level);
        C[8] = Process1Level(C[8], level);

        S[9] = (C2[6] + C2[7] + C2[8]);
		C[9] = C_AND(C2[6], C2[7], level)
              + C_AND(C2[6], C2[8], level) + C_AND(C2[7], C2[8], level);
        C[9] = Process1Level(C[9], level);

        S[10] = (C2[9] + C2[10]);
		C[10] = C_AND(C2[9], C2[10], level);
        C[10] = Process1Level(C[10], level);

        /////////////////////////////////////

        C[11] = C_AND(S[0], S[1], level)
              + C_AND(S[0], S[2], level) + C_AND(S[1], S[2], level);
        C[11] = Process1Level(C[11], level);
        S[0] = (S[0] + S[1] + S[2]);

        C[12] = C_AND(S[3], S[4], level)
              + C_AND(S[3], S[5], level) + C_AND(S[4], S[5], level);
        C[12] = Process1Level(C[12], level);
        S[1] = (S[3] + S[4] + S[5]);

        C[13] = C_AND(S[6], S[7], level)
              + C_AND(S[6], S[8], level) + C_AND(S[7], S[8], level);
        C[13] = Process1Level(C[13], level);
        S[2] = (S[6] + S[7] + S[8]);

        C[14] = C_AND(S[9], S[10], level);
        C[14] = Process1Level(C[14], level);
        S[3] = (S[9] + S[10]);

        ////////////////////////////////////

        C[15] = C_AND(S[0], S[1], level);
        C[15] = Process1Level(C[15], level);
        S[0] = (S[0] + S[1]);

        C[16] = C_AND(S[2], S[3], level);
        C[16] = Process1Level(C[16], level);
        S[1] = (S[2] + S[3]);

        ////////////////////////////////////

        C[17] = C_AND(S[0], S[1], level);
        C[17] = Process1Level(C[17], level);
        result.b[1] = (S[0] + S[1]);

        level++;


        // Depth 3
		for(int i=0; i<18; i+=3)
		{
            int j = i/3;
			S[j] = (C[i] + C[i+1] + C[i+2]);
			C[j] = C_AND(C[i], C[i+1], level)
                 + C_AND(C[i], C[i+2], level) + C_AND(C[i+1], C[i+2], level);
            C[j] = Process1Level(C[j], level);
		}

		C[6] = C_AND(S[0], S[1], level)
              + C_AND(S[0], S[2], level) + C_AND(S[1], S[2], level);
        C[6] = Process1Level(C[6], level);
        S[0] = (S[0] + S[1] + S[2]);

        C[7] = C_AND(S[3], S[4], level)
              + C_AND(S[3], S[5], level) + C_AND(S[4], S[5], level);
        C[7] = Process1Level(C[7], level);
        S[1] = (S[3] + S[4] + S[5]);

        C[8] = C_AND(S[0], S[1], level);
        C[8] = Process1Level(C[8], level);
        result.b[2] = (S[0] + S[1]);

        level++;


        // Depth 4
		for(int i=0; i<9; i+=3)
		{
            int j = i/3;
			S[j] = (C[i] + C[i+1] + C[i+2]);
			C[j] = C_AND(C[i], C[i+1], level)
                 + C_AND(C[i], C[i+2], level) + C_AND(C[i+1], C[i+2], level);
            C[j] = Process1Level(C[j], level);
		}

		C[3] = C_AND(S[0], S[1], level)
              + C_AND(S[0], S[2], level) + C_AND(S[1], S[2], level);
        C[3] = Process1Level(C[3], level);
        result.b[3] = (S[0] + S[1] + S[2]);

        level++;

        // Depth 5
        S[0] = (C[0] + C[1]);
        C[0] = C_AND(C[0], C[1], level);
        S[1] = (C[2] + C[3]);
        C[1] = C_AND(C[2], C[3], level);

        result.b[4] = (S[0] + S[1]);

        C[2] = C_AND(S[0], S[1], level);

        result.b[5] = C[0] + C[1] + C[2];
        result.b[5] = Process1Level(result.b[5], level);

        level++;

        Update_ValueLevel(result.b[4], level);
        Update_ValueLevel(result.b[3], level);
        Update_ValueLevel(result.b[2], level);
        Update_ValueLevel(result.b[1], level);
		Update_ValueLevel(result.b[0], level);

        result.lvl = level;
    }

	return result;
}

void Instruct::DirectSort(InstByte *out, InstByte *in, int size){

// Construct M
	InstBit** M = new InstBit*[size];
	for(int i=0; i<size; i++)
		M[i] = new InstBit[size];

	for(int i=0; i<size; i++){
		M[i][i].bit = 0;
		M[i][i].lvl = level_isLessThan_Byte;
		for(int j=i+1; j<size; j++){
			M[j][i] = isLessThan(in[i], in[j]);		// ( x[i] < x[j] )
			M[i][j] = M[j][i] + 1;				// ( x[i] < x[j] )'
		}
	}
	Update_GlobalLevel(level_isLessThan_Byte);
	Update_TableLevel();



	// Construct Sigma Vector by summing columns of M
	InstByte *S = new InstByte[size];
	for(int i=0; i<size; i++)
	{
		S[i] = HammingWeight(M[i], size);
	}
	Update_GlobalLevel_HW(size);
	Update_TableLevel();


	// Check Rankings and construct the output vector
	InstByte temp;
	temp.lvl = globalLevel;
    int temp_level = globalLevel + level_isEqual;

	for(int i=0; i<size; i++)
	{
        out[i].lvl = temp_level;

		temp = InttoInstByte(i);

		for(int j=0; j<size; j++)
		{
			InstBit z = isEqual(temp, S[j]);
			InstByte k = in[j];
			Update_ValueLevel(k, temp_level);
			out[i] = out[i] + C_AND(z, k, temp_level);

		}
		out[i] = Process1Level(out[i], temp_level);

	}
	Update_GlobalLevel(level_isEqual);
	Update_GlobalLevel(level_and);
	Update_TableLevel();
}

void Instruct::DirectSort(InstInt *out, InstInt *in, int size){

// Construct M
	InstBit** M = new InstBit*[size];
	for(int i=0; i<size; i++)
		M[i] = new InstBit[size];

	for(int i=0; i<size; i++){
		M[i][i].bit = 0;
		M[i][i].lvl = level_isLessThan_Int;
		for(int j=i+1; j<size; j++){
			M[j][i] = isLessThan(in[i], in[j]);		// ( x[i] < x[j] )
			M[i][j] = M[j][i] + 1;				// ( x[i] < x[j] )'
		}
	}
	Update_GlobalLevel(level_isLessThan_Int);
	Update_TableLevel();



	// Construct Sigma Vector by summing columns of M
	InstByte *S = new InstByte[size];
	for(int i=0; i<size; i++)
	{
		S[i] = HammingWeight(M[i], size);
	}
	Update_GlobalLevel_HW(size);
	Update_TableLevel();


	// Check Rankings and construct the output vector
	InstByte temp;
	temp.lvl = globalLevel;
    int temp_level = globalLevel + level_isEqual;

	for(int i=0; i<size; i++)
	{
        out[i].lvl = temp_level;

		temp = InttoInstByte(i);

		for(int j=0; j<size; j++)
		{
			InstBit z = isEqual(temp, S[j]);
			InstInt k = in[j];
			Update_ValueLevel(k, temp_level);
			out[i] = out[i] + C_AND(z, k, temp_level);

		}
		out[i] = Process1Level(out[i], temp_level);

	}
	Update_GlobalLevel(level_isEqual);
	Update_GlobalLevel(level_and);
	Update_TableLevel();
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
void Instruct::BatcherSort(InstByte *out, InstByte *in, int size){

    int lvl = globalLevel;
    int N = size;

    int level = 1;
    for (int p = 1; p < N; p += p)                          // p = {1, 2, 4, ... }
    {
        for (int k = p; k > 0; k /= 2)                      // k = {{1}, {2, 1}, {4, 2, 1}, ... }
        {
            lvl = globalLevel;
            for (int j = k%p; j+k < 2*p; j += (k+k))
            {
                if(p>2 && k==p)
                {
                    for (int i = 1; i < k-1; i++)
                    {
                        for(int q = 0; q < N/(p+p); q++)
                        {
                            int index1 = j+i+(2*p*q);
                            int index2 = j+i+(2*p*q)+k;
                            InstBit A = isLessThan(in[index1], in[index2]);
                            InstBit B = A + 1;
                            int temp_lvl = lvl+level_isLessThan_Byte;
                            Update_ValueLevel(in[index1], temp_lvl);
                            Update_ValueLevel(in[index2], temp_lvl);
                            InstByte min = C_AND(A, in[index1], temp_lvl) + C_AND(B, in[index2], temp_lvl);
                            InstByte max = C_AND(A, in[index2], temp_lvl) + C_AND(B, in[index1], temp_lvl);
                            min = Process1Level(min, temp_lvl);
                            max = Process1Level(max, temp_lvl);
                            temp_lvl = temp_lvl+level_and;
                            in[index1] = min;
                            in[index2] = max;
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < k; i++)
                    {
                        for(int q = 0; q < N/(p+p); q++)
                        {
                            int index1 = j+i+(2*p*q);
                            int index2 = j+i+(2*p*q)+k;
                            InstBit A = isLessThan(in[index1], in[index2]);
                            InstBit B = A + 1;
                            int temp_lvl = lvl+level_isLessThan_Byte;
                            Update_ValueLevel(in[index1], temp_lvl);
                            Update_ValueLevel(in[index2], temp_lvl);
                            InstByte min = C_AND(A, in[index1], temp_lvl) + C_AND(B, in[index2], temp_lvl);
                            InstByte max = C_AND(A, in[index2], temp_lvl) + C_AND(B, in[index1], temp_lvl);
                            min = Process1Level(min, temp_lvl);
                            max = Process1Level(max, temp_lvl);
                            temp_lvl = temp_lvl+level_and;
                            in[index1] = min;
                            in[index2] = max;
                        }
                    }
                }
            }
            if(k==p/2 && p!=N/2)
            {
                for(int q = 0; q < N; q+=(4*p))
                {
                    int index1 = q;
                    int index2 = q+2*p;
                    InstBit A = isLessThan(in[index1], in[index2]);
                    InstBit B = A + 1;
                    int temp_lvl = lvl+level_isLessThan_Byte;
                    Update_ValueLevel(in[index1], temp_lvl);
                    Update_ValueLevel(in[index2], temp_lvl);
                    InstByte min = C_AND(A, in[index1], temp_lvl) + C_AND(B, in[index2], temp_lvl);
                    InstByte max = C_AND(A, in[index2], temp_lvl) + C_AND(B, in[index1], temp_lvl);
                    min = Process1Level(min, temp_lvl);
                    max = Process1Level(max, temp_lvl);
                    temp_lvl += level_and;
                    in[index1] = min;
                    in[index2] = max;
                    ////////////////////////////////////////
                    index1 += 2*p-1;
                    index2 += 2*p-1;
                    A = isLessThan(in[index1], in[index2]);
                    B = A + 1;
                    temp_lvl -= level_and;
                    Update_ValueLevel(in[index1], temp_lvl);
                    Update_ValueLevel(in[index2], temp_lvl);
                    min = C_AND(A, in[index1], temp_lvl) + C_AND(B, in[index2], temp_lvl);
                    max = C_AND(A, in[index2], temp_lvl) + C_AND(B, in[index1], temp_lvl);
                    min = Process1Level(min, temp_lvl);
                    max = Process1Level(max, temp_lvl);
                    temp_lvl = temp_lvl+level_and;
                    in[index1] = min;
                    in[index2] = max;
                }
            }
            level++;
            Update_GlobalLevel(level_isLessThan_Byte);
			Update_GlobalLevel(level_and);
			Update_TableLevel();
			for(int i=0; i<(size); i++)
                Update_ValueLevel(in[i], globalLevel);
        }
    }

    for(int i=0; i<(size); i++)
		out[i] = in[i];
}

void Instruct::BatcherSort(InstInt *out, InstInt *in, int size){

    int lvl = globalLevel;
    int N = size;

    int level = 1;
    for (int p = 1; p < N; p += p)                          // p = {1, 2, 4, ... }
    {
        for (int k = p; k > 0; k /= 2)                      // k = {{1}, {2, 1}, {4, 2, 1}, ... }
        {
            lvl = globalLevel;
            for (int j = k%p; j+k < 2*p; j += (k+k))
            {
                if(p>2 && k==p)
                {
                    for (int i = 1; i < k-1; i++)
                    {
                        for(int q = 0; q < N/(p+p); q++)
                        {
                            int index1 = j+i+(2*p*q);
                            int index2 = j+i+(2*p*q)+k;
                            InstBit A = isLessThan(in[index1], in[index2]);
                            InstBit B = A + 1;
                            int temp_lvl = lvl+level_isLessThan_Int;
                            Update_ValueLevel(in[index1], temp_lvl);
                            Update_ValueLevel(in[index2], temp_lvl);
                            InstInt min = C_AND(A, in[index1], temp_lvl) + C_AND(B, in[index2], temp_lvl);
                            InstInt max = C_AND(A, in[index2], temp_lvl) + C_AND(B, in[index1], temp_lvl);
                            min = Process1Level(min, temp_lvl);
                            max = Process1Level(max, temp_lvl);
                            temp_lvl = temp_lvl+level_and;
                            in[index1] = min;
                            in[index2] = max;
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < k; i++)
                    {
                        for(int q = 0; q < N/(p+p); q++)
                        {
                            int index1 = j+i+(2*p*q);
                            int index2 = j+i+(2*p*q)+k;
                            InstBit A = isLessThan(in[index1], in[index2]);
                            InstBit B = A + 1;
                            int temp_lvl = lvl+level_isLessThan_Int;
                            Update_ValueLevel(in[index1], temp_lvl);
                            Update_ValueLevel(in[index2], temp_lvl);
                            InstInt min = C_AND(A, in[index1], temp_lvl) + C_AND(B, in[index2], temp_lvl);
                            InstInt max = C_AND(A, in[index2], temp_lvl) + C_AND(B, in[index1], temp_lvl);
                            min = Process1Level(min, temp_lvl);
                            max = Process1Level(max, temp_lvl);
                            temp_lvl = temp_lvl+level_and;
                            in[index1] = min;
                            in[index2] = max;
                        }
                    }
                }
            }
            if(k==p/2 && p!=N/2)
            {
                for(int q = 0; q < N; q+=(4*p))
                {
                    int index1 = q;
                    int index2 = q+2*p;
                    InstBit A = isLessThan(in[index1], in[index2]);
                    InstBit B = A + 1;
                    int temp_lvl = lvl+level_isLessThan_Int;
                    Update_ValueLevel(in[index1], temp_lvl);
                    Update_ValueLevel(in[index2], temp_lvl);
                    InstInt min = C_AND(A, in[index1], temp_lvl) + C_AND(B, in[index2], temp_lvl);
                    InstInt max = C_AND(A, in[index2], temp_lvl) + C_AND(B, in[index1], temp_lvl);
                    min = Process1Level(min, temp_lvl);
                    max = Process1Level(max, temp_lvl);
                    temp_lvl += level_and;
                    in[index1] = min;
                    in[index2] = max;
                    ////////////////////////////////////////
                    index1 += 2*p-1;
                    index2 += 2*p-1;
                    A = isLessThan(in[index1], in[index2]);
                    B = A + 1;
                    temp_lvl -= level_and;
                    Update_ValueLevel(in[index1], temp_lvl);
                    Update_ValueLevel(in[index2], temp_lvl);
                    min = C_AND(A, in[index1], temp_lvl) + C_AND(B, in[index2], temp_lvl);
                    max = C_AND(A, in[index2], temp_lvl) + C_AND(B, in[index1], temp_lvl);
                    min = Process1Level(min, temp_lvl);
                    max = Process1Level(max, temp_lvl);
                    temp_lvl = temp_lvl+level_and;
                    in[index1] = min;
                    in[index2] = max;
                }
            }
            level++;
            Update_GlobalLevel(level_isLessThan_Int);
			Update_GlobalLevel(level_and);
			Update_TableLevel();
			for(int i=0; i<(size); i++)
                Update_ValueLevel(in[i], globalLevel);
        }
    }

    for(int i=0; i<(size); i++)
		out[i] = in[i];
}

void Instruct::BubbleSort(InstByte *out, InstByte *in, int size){
	if(size == 0 || size == 1)
		return;

	int lvl;
	bool loop = false;
	// Size = N passes over the array

	for(int k=0; k<size/2; k++){
		// Depth+=1
		// Starting from the first element : even lines of the figure
		lvl = globalLevel;
		loop = false;
		for(int i=0; i<(size-1); i+=2){
			//unsigned int z, min, max;

			// Compare
			//LessThan(z, A[i], A[i+1]);

			InstBit A = isLessThan(in[i], in[i+1]);
			InstBit B = A + 1;

			// B[i] = (A[i]<A[i+1]).A[i] + (A[i]<A[i+1])'.A[i+1]
			//min = C_AND(z, A[i]) + C_AND((z+1)%2, A[i+1]);
			// B[i+1] = (A[i]<A[i+1]).A[i+1] + (A[i]<A[i+1])'.A[i]
			//max = C_AND(z, A[i+1]) + C_AND((z+1)%2, A[i]);

			int temp_lvl = lvl+level_isLessThan_Byte;
			Update_ValueLevel(in[i], temp_lvl);
			Update_ValueLevel(in[i+1], temp_lvl);

			//////////////////////////////////

			//out[0]
			InstByte min = C_AND(A, in[i], temp_lvl) + C_AND(B, in[i+1], temp_lvl);
			//out[1]
			InstByte max = C_AND(A, in[i+1], temp_lvl) + C_AND(B, in[i], temp_lvl);


			min = Process1Level(min, temp_lvl);
			max = Process1Level(max, temp_lvl);
			temp_lvl = temp_lvl+level_and;

			// Swap
			//A[i] = min;
			//A[i+1] = max;
			in[i] = min;
			in[i+1] = max;
			loop = true;
		}

		if(loop == true){
			Update_GlobalLevel(level_isLessThan_Byte);
			Update_GlobalLevel(level_and);
			Update_TableLevel();
		}
		for(int i=0; i<(size); i++)
			Update_ValueLevel(in[i], globalLevel);

		// Depth+=1
		// Starting from the second element : odd lines of the figure
		loop = false;
		lvl = globalLevel;
		for(int i=1; i<(size-2); i+=2){
			//unsigned int z, min, max;

			// Compare
			//LessThan(z, A[i], A[i+1]);

			InstBit A = isLessThan(in[i], in[i+1]);
			InstBit B = A + 1;

			// B[i] = (A[i]<A[i+1]).A[i] + (A[i]<A[i+1])'.A[i+1]
			//min = C_AND(z, A[i]) + C_AND((z+1)%2, A[i+1]);
			// B[i+1] = (A[i]<A[i+1]).A[i+1] + (A[i]<A[i+1])'.A[i]
			//max = C_AND(z, A[i+1]) + C_AND((z+1)%2, A[i]);

			int temp_lvl = lvl+level_isLessThan_Byte;
			Update_ValueLevel(in[i], temp_lvl);
			Update_ValueLevel(in[i+1], temp_lvl);

			//////////////////////////////////

			//out[0]
			InstByte min = C_AND(A, in[i], temp_lvl) + C_AND(B, in[i+1], temp_lvl);
			//out[1]
			InstByte max = C_AND(A, in[i+1], temp_lvl) + C_AND(B, in[i], temp_lvl);


			min = Process1Level(min, temp_lvl);
			max = Process1Level(max, temp_lvl);
			temp_lvl = temp_lvl+level_and;


			in[i] = min;
			in[i+1] = max;
			loop = true;
		}

		if(loop == true){
			Update_GlobalLevel(level_isLessThan_Byte);
			Update_GlobalLevel(level_and);
			Update_TableLevel();
		}

		for(int i=0; i<(size); i++)
			Update_ValueLevel(in[i], globalLevel);
	}

	for(int i=0; i<(size); i++)
		out[i] = in[i];
}

void Instruct::BubbleSort(InstInt *out, InstInt *in, int size){
	if(size == 0 || size == 1)
		return;

	int lvl;
	bool loop = false;
	// Size = N passes over the array

	for(int k=0; k<size/2; k++){
		// Depth+=1
		// Starting from the first element : even lines of the figure
		lvl = globalLevel;
		loop = false;
		for(int i=0; i<(size-1); i+=2){
			//unsigned int z, min, max;

			// Compare
			InstBit A = isLessThan(in[i], in[i+1]);
			InstBit B = A + 1;

			int temp_lvl = lvl+level_isLessThan_Int;
			Update_ValueLevel(in[i], temp_lvl);
			Update_ValueLevel(in[i+1], temp_lvl);

			//////////////////////////////////

			//out[0]
			InstInt min = C_AND(A, in[i], temp_lvl) + C_AND(B, in[i+1], temp_lvl);
			//out[1]
			InstInt max = C_AND(A, in[i+1], temp_lvl) + C_AND(B, in[i], temp_lvl);

			min = Process1Level(min, temp_lvl);
			max = Process1Level(max, temp_lvl);
			temp_lvl = temp_lvl+level_and;

			// Swap
			//A[i] = min;
			//A[i+1] = max;
			in[i] = min;
			in[i+1] = max;
			loop = true;
		}

		if(loop == true){
			Update_GlobalLevel(level_isLessThan_Int);
			Update_GlobalLevel(level_and);
			Update_TableLevel();
		}
		for(int i=0; i<(size); i++)
			Update_ValueLevel(in[i], globalLevel);

		// Depth+=1
		// Starting from the second element : odd lines of the figure
		loop = false;
		lvl = globalLevel;
		for(int i=1; i<(size-2); i+=2){

			// Compare
			InstBit A = isLessThan(in[i], in[i+1]);
			InstBit B = A + 1;

			int temp_lvl = lvl+level_isLessThan_Int;
			Update_ValueLevel(in[i], temp_lvl);
			Update_ValueLevel(in[i+1], temp_lvl);

			//////////////////////////////////
			InstInt min = C_AND(A, in[i], temp_lvl) + C_AND(B, in[i+1], temp_lvl);
			InstInt max = C_AND(A, in[i+1], temp_lvl) + C_AND(B, in[i], temp_lvl);
			min = Process1Level(min, temp_lvl);
			max = Process1Level(max, temp_lvl);
			temp_lvl = temp_lvl+level_and;

			in[i] = min;
			in[i+1] = max;
			loop = true;
		}

		if(loop == true){
			Update_GlobalLevel(level_isLessThan_Int);
			Update_GlobalLevel(level_and);
			Update_TableLevel();
		}

		for(int i=0; i<(size); i++)
			Update_ValueLevel(in[i], globalLevel);
	}

	for(int i=0; i<(size); i++)
		out[i] = in[i];
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

InstBit Instruct::Encrypt(InstBit &x){
	InstBit r;

	r.bit = n->Prim_Encrypt(x.bit, 0);
	r.lvl = 0;

	return r;
}

InstBit Instruct::Decrypt(InstBit &x){
	InstBit r;

	r.bit = n->Prim_Decrypt(x.bit, globalLevel);
	r.lvl = x.lvl;

	return r;
}

InstBit Instruct::Decrypt(InstBit &x, int level){
	InstBit r;

	r.bit = n->Prim_Decrypt(x.bit, level);
	r.lvl = x.lvl;

	return r;
}

InstByte Instruct::Encrypt(InstByte &x){
	InstByte r;

	for(int i=0; i<8; i++){
		r.b[i].bit = n->Prim_Encrypt(x.b[i].bit, 0);
		r.b[i].lvl = 0;
	}

	return r;
}

InstByte Instruct::Decrypt(InstByte &x){
	InstByte r;

	for(int i=0; i<8; i++){
		r.b[i].bit = n->Prim_Decrypt(x.b[i].bit, globalLevel);
		r.b[i].lvl = x.b[i].lvl;
	}

	return r;
}


InstByte Instruct::Encrypt(InstByte &x, int lvl){
	InstByte r;

	for(int i=0; i<8; i++){
		r.b[i].bit = n->Prim_Encrypt(x.b[i].bit, lvl);
		r.b[i].lvl = 0;
	}

	return r;
}

InstByte Instruct::Decrypt(InstByte &x, int lvl){
	InstByte r;

	for(int i=0; i<8; i++){
		r.b[i].bit = n->Prim_Decrypt(x.b[i].bit, lvl);
		r.b[i].lvl = x.b[i].lvl;
	}

	return r;
}

InstInt Instruct::Encrypt(InstInt &x){
	InstInt r;

	for(int i=0; i<32; i++){
		r.b[i].bit = n->Prim_Encrypt(x.b[i].bit, 0);
		r.b[i].lvl = 0;
	}

	return r;
}

InstInt Instruct::Decrypt(InstInt &x){
	InstInt r;

	for(int i=0; i<32; i++){
		r.b[i].bit = n->Prim_Decrypt(x.b[i].bit, globalLevel);
		r.b[i].lvl = x.b[i].lvl;
	}

	return r;
}

InstInt Instruct::Decrypt(InstInt &x, int lvl){
	InstInt r;

	for(int i=0; i<32; i++){
		r.b[i].bit = n->Prim_Decrypt(x.b[i].bit, lvl);
		r.b[i].lvl = x.b[i].lvl;
	}

	return r;
}
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void Instruct::Update_TableLevel(){
	n->Prim_ReduceKeysLevel(globalLevel);
}

void Instruct::Update_ValueLevel(InstBit &x){
	int l = globalLevel - x.lvl;
	for(int i=0; i<l; i++)
		x.bit = n->Prim_ModSwitch(x.bit, x.lvl+i);
	x.lvl = globalLevel;
}

void Instruct::Update_ValueLevel(InstByte &x){
	int xl = x.b[0].lvl;
	int l = globalLevel - xl;
	for(int i=0; i<l; i++)
		for(int j=0; j<8; j++)
			x.b[j].bit = n->Prim_ModSwitch(x.b[j].bit, xl+i);

	for(int j=0; j<8; j++)
		x.b[j].lvl = globalLevel;

	x.lvl = globalLevel;
}

void Instruct::Update_ValueLevel(InstInt &x){
	int xl = x.b[0].lvl;
	int l = globalLevel - xl;
	for(int i=0; i<l; i++)
		for(int j=0; j<32; j++)
			x.b[j].bit = n->Prim_ModSwitch(x.b[j].bit, xl+i);

	for(int j=0; j<32; j++)
		x.b[j].lvl = globalLevel;

	x.lvl = globalLevel;
}

void Instruct::Update_ValueLevel(InstBit &x, int lvl){
	int l = lvl - x.lvl;
	if(l != 0)
		n->Arith_CoeffReduce(x.bit, x.bit, x.lvl);
	for(int i=0; i<l; i++){
		x.bit = n->Prim_ModSwitch_LongPoly(x.bit, deg(x.bit)+1,x.lvl+i);
	}
	x.lvl = lvl;
}

void Instruct::Update_ValueLevel(InstByte &x, int lvl){
	int xl = x.b[0].lvl;
	int l = lvl - xl;
	for(int i=0; i<l; i++)
		for(int j=0; j<8; j++)
			x.b[j].bit = n->Prim_ModSwitch_LongPoly(x.b[j].bit, deg(x.b[j].bit)+1,xl+i);

	for(int j=0; j<8; j++)
		x.b[j].lvl = lvl;

	x.lvl = lvl;
}

void Instruct::Update_ValueLevel(InstInt &x, int lvl){
	int xl = x.b[0].lvl;
	int l = lvl - xl;
	for(int i=0; i<l; i++)
		for(int j=0; j<32; j++)
			x.b[j].bit = n->Prim_ModSwitch_LongPoly(x.b[j].bit, deg(x.b[j].bit)+1,xl+i);

	for(int j=0; j<32; j++)
		x.b[j].lvl = lvl;

	x.lvl = lvl;
}

void Instruct::Update_GlobalLevel(InstLevel a){
//	if(a == level_isEqual)
//		globalLevel = globalLevel + level_isEqual;
//	else if(a == level_isLessThan_Byte)
//			globalLevel = globalLevel + level_isLessThan_Byte;
//	else if(a == level_isLessThan_Int)
//		globalLevel = globalLevel + level_isLessThan_Int;
//	else if(a == level_Add)
//			globalLevel = globalLevel + level_Add;
//	else if(a == level_and)
//			globalLevel = globalLevel + level_and;

    globalLevel = globalLevel + a;
}

void Instruct::Update_GlobalLevel_HW(int size){
	if(size == 4)
		globalLevel = globalLevel + level_HW_4;
	else if(size == 8)
        globalLevel = globalLevel + level_HW_8;
    else if(size == 16)
        globalLevel = globalLevel + level_HW_16;
    else if(size == 32)
        globalLevel = globalLevel + level_HW_32;
    else if(size == 64)
        globalLevel = globalLevel + level_HW_64;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

InstBit Instruct::Process1Level(InstBit &x, int lvl){
	InstBit r;
	r.bit = x.bit;
	n->Arith_PolyCoeffReduce(r.bit, r.bit, lvl);
	r.bit = n->Prim_RelinRingFFT(r.bit, lvl);
	n->Arith_PolyCoeffReduce(r.bit, r.bit, lvl);
	r.bit = n->Prim_ModSwitch(r.bit, lvl);
	r.lvl = x.lvl + 1;
	return r;
}

InstByte Instruct::Process1Level(InstByte &x, int lvl){
	InstByte r;

	for(int i=0; i<8; i++)
		r.b[i] = Process1Level(x.b[i], lvl);

	r.lvl = x.lvl + 1;
	return r;
}

InstInt Instruct::Process1Level(InstInt &x, int lvl){
	InstInt r;

	for(int i=0; i<32; i++)
		r.b[i] = Process1Level(x.b[i], lvl);

	r.lvl = x.lvl + 1;
	return r;
}

void  Instruct::Reduce(InstBit &x, int level){
	n->Arith_PolyCoeffReduce(x.bit, x.bit, level);
}

int Instruct::EqualizeLevel(InstBit &x, InstBit &y){

	int lvl_x = x.lvl;
	int lvl_y = y.lvl;

	if(lvl_x == lvl_y)
		return lvl_x;
	else if(lvl_x > lvl_y){

		int t = lvl_x - lvl_y;
		for(int i=0; i<t; i++)
			y.bit = n->Prim_ModSwitch(y.bit, lvl_y+i);

		y.lvl = y.lvl + t;
		return lvl_x;
	}
	else{
		int t = lvl_y - lvl_x;
		for(int i=0; i<t; i++)
			x.bit = n->Prim_ModSwitch(x.bit, lvl_x+i);

		x.lvl = x.lvl + t;
		return lvl_y;
	}
}

int Instruct::EqualizeLevel(InstBit &x, InstInt &y){

	int lvl_x = x.lvl;
	int lvl_y = y.lvl;

	if(lvl_x == lvl_y)
		return lvl_x;
	else if(lvl_x > lvl_y){
        int t = lvl_x - lvl_y;
		for(int i=0; i<t; i++)
			for(int j=0; j<32; j++)
				y.b[j].bit = n->Prim_ModSwitch(y.b[j].bit, lvl_y+i);

		for(int j=0; j<32; j++)
			y.b[j].lvl = y.b[j].lvl + t;

		return lvl_x;
	}
	else{
		int t = lvl_y - lvl_x;
		for(int i=0; i<t; i++)
			x.bit = n->Prim_ModSwitch(x.bit, lvl_x+i);

		x.lvl = x.lvl + t;
		return lvl_y;
	}
}


int Instruct::EqualizeLevel(InstByte &x, InstByte &y){

	int lvl_x = x.b[0].lvl;
	int lvl_y = y.b[0].lvl;

	if(lvl_x == lvl_y)
		return lvl_x;
	else if(lvl_x > lvl_y){
		int t = lvl_x - lvl_y;
		for(int i=0; i<t; i++)
			for(int j=0; j<8; j++)
				y.b[j].bit = n->Prim_ModSwitch(y.b[j].bit, lvl_y+i);

		for(int j=0; j<8; j++)
			y.b[j].lvl = y.b[j].lvl + t;
		return lvl_x;
	}
	else{
		int t = lvl_y - lvl_x;
		for(int i=0; i<t; i++)
			for(int j=0; j<8; j++)
				x.b[j].bit = n->Prim_ModSwitch(x.b[j].bit, lvl_x+i);

		for(int j=0; j<8; j++)
			x.b[j].lvl = x.b[j].lvl + t;
		return lvl_y;
	}
}

int Instruct::EqualizeLevel(InstInt &x, InstInt &y){

	int lvl_x = x.b[0].lvl;
	int lvl_y = y.b[0].lvl;

	if(lvl_x == lvl_y)
		return lvl_x;
	else if(lvl_x > lvl_y){

		int t = lvl_x - lvl_y;
		for(int i=0; i<t; i++)
			for(int j=0; j<32; j++)
				y.b[j].bit = n->Prim_ModSwitch(y.b[j].bit, lvl_y+i);

		for(int j=0; j<32; j++)
			y.b[j].lvl = y.b[j].lvl + t;

		return lvl_x;
	}
	else{
		int t = lvl_y - lvl_x;
		for(int i=0; i<t; i++)
			for(int j=0; j<32; j++)
				x.b[j].bit = n->Prim_ModSwitch(x.b[j].bit, lvl_x+i);

		for(int j=0; j<32; j++)
			x.b[j].lvl = x.b[j].lvl + t;

		return lvl_y;
	}
}
