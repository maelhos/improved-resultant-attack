#ifndef ANEMOI__H
#define ANEMOI__H

#include <NTL/lzz_p.h>
#include <openssl/sha.h>

NTL_CLIENT

#define bSEED "ANEMOI"
#define ALPHA 3

class anemoi
{
public:
    long nr;
    zz_p a;
    zz_p a_inv;
    Vec<zz_p> ci;
    Vec<zz_p> di;

    zz_p x;
    zz_p y;

    long beta;

    anemoi(long nr);

    void R(zz_p& x, zz_p& y, long r);
    void H(zz_p& x, zz_p& y);
    void M(zz_p& x, zz_p& y);
    void encrypt(zz_p& x, zz_p& y);

    void R_inv(zz_p& x, zz_p& y, long r);
    void H_inv(zz_p& x, zz_p& y);
    void M_inv(zz_p& x, zz_p& y);
    void decrypt(zz_p& x, zz_p& y);

    // special just for the attack
    void decrypt_noM_r(zz_p& ox, zz_p& oy, const zz_p x, const zz_p y, long r);

    ~anemoi() {}
};

#endif