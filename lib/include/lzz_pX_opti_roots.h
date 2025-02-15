#ifndef LZZ_PX_OPTI_ROOTS__H
#define LZZ_PX_OPTI_ROOTS__H

#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pXFactoring.h>
#include <iomanip>
#define get_time GetWallTime
NTL_CLIENT

vec_zz_p zz_pX_roots_opti(const zz_pX& P);

#endif