#ifndef CIPHERTEXT_H_
#define CIPHERTEXT_H_

#include "general.h"

class CipherText
{
public:
            CipherText(){set_d(to_ZZ(1)); set_ct(0); set_fresh(true); set_level(0); set_multcnt(0);}
            // Operator Overloading
            CipherText& operator=(const CipherText &rhs);
            CipherText& operator+=(const CipherText &rhs);
            const CipherText operator+(const CipherText& rhs) const;
            CipherText& operator-=(const CipherText &rhs);
            const CipherText operator-(const CipherText& rhs) const;

            CipherText& operator+=(const ZZ &rhs);
            const CipherText operator+(const ZZ& rhs) const;
            CipherText& operator*=(const ZZ &rhs);
            const CipherText operator*(const ZZ& rhs) const;

            friend std::ostream& operator<< (std::ostream&, const CipherText& ct);

            const ZZX&  ct() const{ return ct_; };
            const ZZ&   d() const{ return d_; };
            bool        fresh()const{return fresh_;};
            int         level()const{return level_;};
            int         multcnt()const{return multcnt_;};

            void    set_ct(int i){ct_ = i;};
            void    set_ct(const ZZX &ct){ct_ = ct;};
            void    set_fresh(bool f){fresh_ = f;};
            void    set_level(int l){level_ = l;};
            void    set_multcnt(int cnt){multcnt_ = cnt;};
            void    set_d(const ZZ &d){d_ = d;};

            void    coeff_reduce(const ZZ &mod);
            //void    poly_reduce(RingType type, const ZZX &mod);
            //void    poly_reduce(ReductionType type, const ZZX &mod, const vec_ZZX &helpers);

private:
            ZZX ct_;
            bool fresh_;
            int level_;
            int multcnt_;
            ZZ  d_;
};

#endif
