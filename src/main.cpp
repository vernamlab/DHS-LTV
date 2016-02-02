#include "homlib.h"

using namespace std;
using namespace NTL;



int main()
{
    Crypto C;

    ZZX m, b;

    CipherText c, result;

    m = 1;
    C.Encrypt(c, m);

    C.Mult(c, c);

    C.Decrypt(b, c);
    cout << b << endl;


    return 0;
}
