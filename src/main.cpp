#include "homlib.h"

using namespace std;
using namespace NTL;



int main()
{
    LaSH L;

    ZZX m, b;

    CipherText c, result;

    m = 1;
    L.Encrypt(c, m);

    L.Mult(c, c);

    L.Decrypt(b, c);
    cout << b << endl;


    return 0;
}
