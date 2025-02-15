#include "anemoi.h"

anemoi::anemoi(long nr) : nr(nr), ci{}, di{} 
{
    unsigned long hash_buff[4];
    ci.SetLength(nr);
    di.SetLength(nr);

    string tmp_str = bSEED;
    SHA256((const unsigned char*)tmp_str.c_str(), tmp_str.length(), (unsigned char*)hash_buff);
    a = zz_p(hash_buff[0] % zz_p::modulus());

    for (long r = 0; r < nr; r++)
    {
        tmp_str = bSEED + to_string(r);
        SHA256((const unsigned char*)tmp_str.c_str(), tmp_str.length(), (unsigned char*)hash_buff);
        ci[r] = zz_p(hash_buff[0] % zz_p::modulus());

        tmp_str = to_string(r) + bSEED;
        SHA256((const unsigned char*)tmp_str.c_str(), tmp_str.length(), (unsigned char*)hash_buff);
        di[r] = zz_p(hash_buff[0] % zz_p::modulus());
    }
    a_inv = zz_p(1) / a;

    zz_pPush push(zz_p::modulus() - 1);
    beta = (zz_p(1) / zz_p(ALPHA))._zz_p__rep;
};

void anemoi::R(zz_p& x, zz_p& y, long r){
    x += ci[r];
    y += di[r];
    M(x, y);
    H(x, y);
}

void anemoi::H(zz_p& x, zz_p& y){
    x -= a*y*y + a_inv;
    y -= power(x, beta);
    x += a*y*y;
}

void anemoi::M(zz_p& x, zz_p& y){
    y += x;
    x += y;
}

void anemoi::encrypt(zz_p& x, zz_p& y){
    for (long i = 0; i < nr; i++)
        R(x, y, i);
    M(x, y);
}

// decrypt 

void anemoi::R_inv(zz_p& x, zz_p& y, long r){
    H_inv(x, y);
    M_inv(x, y);
    x -= ci[r];
    y -= di[r];
}

void anemoi::H_inv(zz_p& x, zz_p& y){
    x -= a*y*y;
    y += power(x, beta);
    x += a*y*y + a_inv;
}

void anemoi::M_inv(zz_p& x, zz_p& y){
    x -= y;
    y -= x;
}

void anemoi::decrypt(zz_p& x, zz_p& y){
    M_inv(x, y);
    for (long i = nr - 1; i > -1; i--)
        R_inv(x, y, i);
}
void anemoi::decrypt_noM_r(zz_p& ox, zz_p& oy, const zz_p x, const zz_p y, long r){
    ox = x;
    oy = y;
    for (long i = r - 1; i > -1; i--)
        R_inv(ox, oy, i);
}