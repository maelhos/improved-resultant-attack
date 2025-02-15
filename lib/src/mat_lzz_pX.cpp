#include "mat_lzz_pX.h"

void mul_aux(Vec<zz_pX>& x, const Mat<zz_p>& A, const Vec<zz_pX>& b)  
{  
    if (&x == &b)
    {
        Vec<zz_pX> xs{};
        mul_aux(xs, A, b);
        x = xs;
        return;
    }
    
    long n = A.NumRows();  
    long l = A.NumCols();  
    
    if (l != b.length())  
        LogicError("matrix mul: dimension mismatch");  
    
    x.SetLength(n);  
    
    long i, k;  
    zz_pX acc, tmp;  
    
    for (i = 1; i <= n; i++) {  
        acc.zero(); 
        for (k = 1; k <= l; k++) {  
            mul(tmp, b(k), A(i, k));  
            add(acc, acc, tmp);  
        }  
        x(i) = acc;  
    }  
}  
  