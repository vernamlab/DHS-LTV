#include "def.h"
#include "crt.h"
#include "ntru.h"

#include <string>
#include <sstream>
#include <fstream>

#ifndef INSTRUCT_H_
#define INSTRUCT_H_


enum InstLevel {level_and = 1, level_isEqual = 3, level_isLessThan_Byte = 4,
level_isLessThan_Int = 6, level_Add = 4, level_HW_4 = 1, level_HW_8 = 2, level_HW_16 = 3,
level_HW_32 = 4, level_HW_64 = 5};

//Instrucion Bit
struct InstBit{
	ZZX bit;
	int lvl;
	InstBit(){
        lvl = 0;
		bit = 0;
	}
	struct InstBit& operator=(const InstBit& other){
		bit = other.bit;
		lvl = other.lvl;
		return *this;
	}

	struct InstBit operator+(const InstBit& rhs){
		InstBit r;
		r.bit = bit + rhs.bit;
		r.lvl = lvl;
		return r;
	}

	struct InstBit operator+(const int& rhs){
			InstBit r;
			r.bit = bit + rhs;
			r.lvl = lvl;
			return r;
	}
};

//Instrucion Byte
struct InstByte{
	InstBit b[8];
	int lvl;
	InstByte(){
        lvl = 0;
	//	for(int i=0; i<8; i++)
	//		b[i].bit = 0;
	}
	struct InstByte operator+(const InstByte& rhs){
		InstByte r;

		for(int i=0; i<8; i++)
			r.b[i] = b[i] + rhs.b[i];

		r.lvl = lvl;
		return r;
	}

};

//Instrucion Int
struct InstInt{
	InstBit b[32];
	int lvl;
	InstInt(){
        lvl = 0;
		//for(int i=0; i<32; i++)
		//	b[i].bit = 0;
	}
	struct InstInt operator+(const InstInt& rhs){
		InstInt r;

		for(int i=0; i<32; i++)
			r.b[i] = b[i] + rhs.b[i];

		r.lvl = lvl;
		return r;
	}
};


class Instruct{
public:
	Instruct(int cyc_deg, int num_primes, int max_prime);
	void GenerateParameters();
	void SaveParameters(string in);
	void LoadParameters(string in);

	void GenerateKeys();
	void SaveKeys(string s);
	void LoadKeys(string s);

	///////////
	// Basic Functions
	InstBit 	Encrypt(InstBit &x);
	InstBit 	Decrypt(InstBit &x);

	InstByte 	Encrypt(InstByte &x);
	InstByte	Decrypt(InstByte &x);

	InstInt 	Encrypt(InstInt &x);
	InstInt 	Decrypt(InstInt &x);
	InstInt 	Decrypt(InstInt &x, int lvl);

	InstByte 	Encrypt(InstByte &x, int lvl);
	InstByte 	Decrypt(InstByte &x, int lvl);

	InstBit 	Decrypt(InstBit &x, int level);

	void	Update_TableLevel();				// Update Table to Global Level
	void	Update_ValueLevel(InstBit &x);		// Update level of the value to Global
	void	Update_ValueLevel(InstByte &x);		// Update level of the value to Global
	void	Update_ValueLevel(InstInt &x);		// Update level of the value to Global

	void	Update_ValueLevel(InstBit &x ,int lvl);		// Update level of the value to given lvl
	void	Update_ValueLevel(InstByte &x,int lvl);		// Update level of the value to given lvl
	void	Update_ValueLevel(InstInt &x,int lvl);		// Update level of the value to given lvl

	void 	Update_GlobalLevel(InstLevel a);	// Update Global Level
	void    Update_GlobalLevel_HW(int size);    // Update Global Level for HW instruction, depending on size
	///////////
	InstBit  	C_AND(InstBit &a, InstBit &b, int level);
	InstByte	C_AND(InstBit &a, InstByte &b, int level);
	InstInt		C_AND(InstBit &a, InstInt &b, int level);

	InstBit  	C_XOR(InstBit &a, InstBit &b);
	InstByte	C_XOR(InstByte &a, InstByte &b);
	InstInt		C_XOR(InstInt &a, InstInt &b);


	InstByte 	InttoInstByte(int i);			// Convert integer to InstrByte
	int			InstBytetoInt(InstByte &i);		// Convert InstrByte to integer

	InstInt 	InttoInstInt(int i);			// Convert integer to InstrByte
	int			InstInttoInt(InstInt &i);		// Convert InstrByte to integer

	// Instructions
	InstBit 	isEqual(InstByte &x, InstByte &y);		//
	InstBit 	isLessThan(InstByte &x, InstByte &y);	//
	InstByte	Add(InstByte &x, InstByte &y);			//

	InstBit 	isLessThan(InstInt &x, InstInt &y);	//
	InstByte	HammingWeight(InstBit *x, int size);	//


	//Algorithms
	void GreedySort(InstByte *out, InstByte *in, int size);
	void GreedySort(InstInt *out, InstInt *in, int size);

	void BubbleSort(InstByte *out, InstByte *in, int size);
	void BubbleSort(InstInt *out, InstInt *in, int size);

	void DirectSort(InstByte *out, InstByte *in, int size);
	void DirectSort(InstInt *out, InstInt *in, int size);

	void BatcherSort(InstByte *out, InstByte *in, int size);
	void BatcherSort(InstInt *out, InstInt *in, int size);

	///////////
	ZZX 		GenRandMess(int size);	// Generate Random Message
	ntru 		*Pointer_NTRU();		// Pointer to Ntru


	InstBit  Process1Level(InstBit &x, int level);
	InstByte Process1Level(InstByte &x, int level);
	InstInt  Process1Level(InstInt &x, int level);

	void  Reduce(InstBit &x, int level);

private:
	ntru 	*n;
	myCRT 	*c;
	int 	globalLevel;

	//PRIVATE FUNC
//	InstBit Process1Level(InstBit &x, int level);		// Do Relinization and swtich to next level
	int 	EqualizeLevel(InstByte &x, InstByte &y);	// Check two values and equalize the level to the larger one
	int 	EqualizeLevel(InstInt &x, InstInt &y);	// Check two values and equalize the level to the larger one
	int 	EqualizeLevel(InstBit &x, InstBit &y);		// Check two values and equalize the level to the larger one
    int 	EqualizeLevel(InstBit &x, InstInt &y);
};


#endif /* INSTRUCT_H_ */




