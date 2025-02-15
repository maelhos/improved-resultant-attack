#ifndef MAT_LZZ_PXI__H
#define MAT_LZZ_PXI__H

#include <NTL/vector.h>
#include "lzz_pXi.h"

NTL_CLIENT

template <auto N>
void mul(Vec<zz_pXi<N>>& x, const Mat<zz_p>& A, const Vec<zz_pXi<N>>& b)  
{  
    if (&x == &b)
    {
        Vec<zz_pXi<N>> xs{};
        mul(xs, A, b);
        x = xs;
        return;
    }
    
    long n = A.NumRows();  
    long l = A.NumCols();  
    
    if (l != b.length())  
        LogicError("matrix mul: dimension mismatch");  
    
    x.SetLength(n);  
    
    long i, k;  
    zz_pXi<N> acc, tmp;  
    
    for (i = 1; i <= n; i++) 
    {  
        acc.zero(); 
        for (k = 1; k <= l; k++) 
        {  
            mul_cst(tmp, b(k), A(i, k));  
            add(acc, acc, tmp);  
        }  
        x(i) = acc;  
    }  
}  
  

#endif