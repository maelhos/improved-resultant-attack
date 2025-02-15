#include "rescue.h"
#include <assert.h>

Mat<zz_p> MDS;
Vec<zz_p> round_constants;

rescue::rescue() : round_constants_counter(0) 
{
    zz_pPush push(zz_p::modulus() - 1);
    inv_d = (zz_p(1) / zz_p(RESCUE_D))._zz_p__rep;
}

void rescue::round(rescue_state& s){
    for (long i = 0; i < RESCUE_T; i++)
        power(s[i], s[i], RESCUE_D);

    mul(s, MDS, s);

    for (long i = 0; i < RESCUE_T; i++)
    {
        s[i] += round_constants[round_constants_counter];
        round_constants_counter++;
    }

    for (long i = 0; i < RESCUE_T; i++)
        power(s[i], s[i], inv_d);

    mul(s, MDS, s);
    
    for (long i = 0; i < RESCUE_T; i++)
    {
        s[i] += round_constants[round_constants_counter];
        round_constants_counter++;
    }
}

void rescue::encrypt(rescue_state& s)
{ 
    for (long r = 0; r < RESCUE_R; r++)
        round(s);
}