#include "fheaes.h"

int main(){
	myTimer t;
	FheAES f;

t.Start();
	f.LTVKeySetUp();
t.Stop();
t.ShowTime("Set Time:\t");

	unsigned char m[16] = {
			0x00, 0x44, 0x88, 0xCC,
			0x11, 0x55, 0x99, 0xDD,
			0x22, 0x66, 0xAA, 0xEE,
			0x33, 0x77, 0xBB, 0xFF
	};
	f.SetMessageBits(m);

	cout << "Encrypt Message Start" << endl;
t.Start();
	f.EncryptMessage();
t.Stop();
t.ShowTime("Enc Mess Time:\t");
	cout << "Encrypt Message Done" << endl;

t.Start();
	f.SetKeys();
t.Stop();
t.ShowTime("Key Set Time:\t");
	cout << "Set Keys Done" << endl;


	f.AESEncryption();
	cout << "AES Encryption Done" << endl;

	return 0;
}
