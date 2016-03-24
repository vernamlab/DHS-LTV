#include "fheprince.h"


//////////////////////////////KEY GENERATION//////////////////////////////////


int main(){

    /////////////////////
    //// Interface   ////
    // A: Plaintext  ////
    // B: Key_0      ////
    // C: Key_1      ////
    /////////////////////
    int A[64],B[64],C[64];
    for (int i=0;i<64;i++){
        A[i] = 0;
        B[i] = 1;
        C[i] = 0;
    }


    FhePRINCE f;
	f.LTVKeySetUp();

    f.SetMessage(A);
    f.SetKeys(B, C);
    f.EncryptMessage();
    f.EncryptKeys();
    f.PRINCEEncryption();
}




