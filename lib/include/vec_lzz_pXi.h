#ifndef VEC_LZZ_PXI__H
#define VEC_LZZ_PXI__H


#include <NTL/vector.h>
#include "lzz_pXi.h"

NTL_CLIENT

template <auto N>
void add(Vec<zz_pXi<N>>& c, const Vec<zz_pXi<N>>& a, const Vec<zz_p>& b)
{
    assert(a.length() == b.length());
    c.SetLength(a.length());
    for (long i = 0; i < a.length(); i++)
        add_cst(c[i], a[i], b[i]);
    
}

#endif