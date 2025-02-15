#ifndef HANKEL_DET__H
#define HANKEL_DET__H
#include "lzz_pXi.h"
#include "lzz_pXi_triangular_ideal.h"
#include <array>

NTL_CLIENT

template <auto alpha, auto N>
void detUis(zz_pXi<N>& det, Vec<zz_pXi<N>>& Uis, zz_pXi_triangular_ideal<alpha, N>& ideal){
    double ti = GetWallTime();
    zz_pXi<N + 1> one{}; one.one(); 

    zz_pXi<N + 1> theta{}; theta.Xn();   
    zz_pXi<N + 1> theta_i{}; theta_i.one();  

    zz_pXi_triangular_ideal<alpha, N + 1> ideal_nth_root = ideal.extend(one);

    zz_pXi<N + 2> det2{};
    det2.one();

    Vec<zz_pXi<N + 1>> coeffs{}; coeffs.SetLength(alpha);
    for (long i = 0; i < alpha; i++)
    {
        coeffs[i].rep.SetLength(i + 1);
        coeffs[i].rep[i] = Uis[alpha - 1 + i];
    }

    long chi_shift;
    Vec<zz_pXi<N + 2>> hs{}; hs.SetLength(alpha);
    for (long j = 0; j < alpha; j++)
    {
        for (long i = 0; i < alpha; i++)
        {
            zz_pXi<N + 2> tmp{};

            chi_shift = (i * j) % alpha;
            tmp.rep.SetLength(chi_shift + 1);
            tmp.rep[chi_shift] = coeffs[i];
            add(hs[j], hs[j], tmp);
        }
    }
    long s = alpha;

    while (s > 1)
    {
        long ns = (s + 1) / 2;
        if (!(s & 1))
            ideal_nth_root.mulReduce(hs[0], hs[0], hs[ns]);

        for (long i = 1; i < ns; i++)
            ideal_nth_root.mulReduce(hs[i], hs[i], hs[s - i]);
        s = ns;
    }

    det = hs[0].rep[0].rep[0] - hs[0].rep[1].rep[0];
    cout << "[+] Det Hankel took : " << GetWallTime() - ti << endl;
}

template <auto N>
void detUis(zz_pXi<N>& det, Vec<zz_pXi<N>>& Uis, zz_pXi_triangular_ideal<3, N>& ideal){
    double ti = GetWallTime();
    zz_pXi<N> a, b, c, d, e;
    a.rep = move(Uis[0].rep);
    b.rep = move(Uis[1].rep);
    c.rep = move(Uis[2].rep);
    d.rep = move(Uis[3].rep);
    e.rep = move(Uis[4].rep);

    zz_pXi<N> t, ad, be{};
    zz_pXi_triangular_ideal<3, N - 1> base = ideal.base;

    base.mulReduce(t, b + a, d + d + e);
    base.mulReduce(ad, a, d);
    base.mulReduce(be, b, e);
    zz_pXi<N> bdae = t - ad - ad - be;
    det = base.mulReduce(c, base.mulReduce(c, c) - bdae) + base.mulReduce(ad, d) + base.mulReduce(b, be);
    cout << "[+] Det Hankel took : " << GetWallTime() - ti << endl;
}

#endif