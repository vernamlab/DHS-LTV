#include "ciphertext.h"

// Operator Overloads
CipherText& CipherText::operator=(const CipherText& rhs)
{
    if(this != &rhs) // If its not a self assignment
    {
        ct_ = rhs.ct();
        level_ = rhs.level();
        fresh_ = rhs.fresh();
        multcnt_ = rhs.multcnt();
    }
    return *this;
}
CipherText& CipherText::operator+=(const CipherText& rhs)
{
    set_ct(ct()+rhs.ct());
    return *this;
}
const CipherText CipherText::operator+(const CipherText& rhs) const
{
    CipherText ct;
    ct = *this;
    ct += rhs;
    return ct;
}
CipherText& CipherText::operator-=(const CipherText& rhs)
{
    set_ct(ct()-rhs.ct());
    return *this;
}
const CipherText CipherText::operator-(const CipherText& rhs) const
{
    CipherText ct;
    ct = *this;
    ct -= rhs;
    return ct;
}
CipherText& CipherText::operator+=(const ZZ& rhs)
{
    set_ct(ct()+rhs);
    return *this;
}
const CipherText CipherText::operator+(const ZZ& rhs) const
{
    CipherText ct;
    ct = *this;
    ct += rhs;
    return ct;
}
CipherText& CipherText::operator*=(const ZZ& rhs)
{
    set_ct(ct()*rhs);
    return *this;
}
const CipherText CipherText::operator*(const ZZ& rhs) const
{
    CipherText ct;
    ct = *this;
    ct *= rhs;
    return ct;
}
std::ostream& operator<<(std::ostream& os, const CipherText& ct)
{
    os  << ct.ct();
    return os;
}

void CipherText::coeff_reduce(const ZZ &mod)
{
    CoeffReduce(ct_, ct_, mod);
}
/*void CipherText::poly_reduce(RingType type, const ZZX &mod)
{
    PolyReduce(ct_, ct_, type, mod);
}*/
//void CipherText::poly_reduce(ReductionType type, const ZZX &mod, const vec_ZZX &helpers)
//{

//}

