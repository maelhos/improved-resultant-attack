#ifndef MAT_LZZ_PXY__H
#define MAT_LZZ_PXY__H

#include <NTL/lzz_pX.h>
NTL_CLIENT

void mul_aux(Vec<zz_pX>& x, const Mat<zz_p>& A, const Vec<zz_pX>& b);

#endif