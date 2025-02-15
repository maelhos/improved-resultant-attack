#ifndef SPECIAL_RESULTANT__H
#define SPECIAL_RESULTANT__H

#include "lzz_pXi.h"
#include "lzz_pXi_triangular_ideal.h"
#include "hankel_det.h"

NTL_CLIENT


template <auto K, auto N> // alpha = K, degree N
void resultant_P_ideal(zz_pXi<N>& r, const zz_pXi<N + 1>& P1, zz_pXi_triangular_ideal<K, N>& ideal){
    long degZ = P1.degXn();

    // for future generalization compatible with alpha > 3
    constexpr long B = K - 1;
    constexpr long D = 2*K - 1;

    Vec<zz_pXi<N>> Uis{};
    long Uis_length = degZ + 1 + 2*B;
    Uis.SetLength(Uis_length);

    for (long i = 0; i <= degZ; i++)
        Uis[i + B] = P1.rep[i];
    
    zz_pXi<N> tmp{};
    for (long i = Uis_length - D - 1; i > -1; i--)
    {
        ideal.base.mulReduce(tmp, ideal.fn, Uis[i + K]);
        add(Uis[i], Uis[i], tmp);
    }
    
    detUis(r, Uis, ideal);
}

#endif