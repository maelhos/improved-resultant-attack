#include <iostream>
#include <NTL/ALL_FEATURES.h>
#include "lzz_pXi_triangular_ideal.h"
#include "anemoi.h"
#include "special_resultant.h"
#include <unistd.h>
#include "lzz_pX_opti_roots.h"

NTL_CLIENT

#define R 8
#define P 0x64ec6dd0392073ULL

template <auto N>
using ideal = zz_pXi_triangular_ideal<ALPHA, N>;

// WARNING, does change X and Y to save memory
template <auto N>
void R_tri(zz_pXi<N+1>& nX, zz_pXi<N+1>& nY, ideal<N>& nI, zz_pXi<N>& X, zz_pXi<N>&Y, ideal<N - 1>& I, const anemoi& A)
{   
    double start;
    // csts
    cout << " - Constants of round " << N - 1 << endl;
    start = get_time();
    add_cst(X, X, A.ci[N - 1]);
    add_cst(Y, Y, A.di[N - 1]);
    cout << "    Time : " << get_time() - start << endl;

    // linear
    cout << " - Linear of round " << N - 1 << endl;
    start = get_time();
    add(Y, Y, X);
    add(X, X, Y);
    cout << "    Time : " << get_time() - start << endl;

    zz_pXi<N> temp{};
    zz_pXi<N+1> Zn{};
    Zn.Xn();

    cout << " - Sbox of round " << N - 1 << endl;
    start = get_time();
    add_cst(X, X, -A.a_inv);
    temp = X - A.a*I.mulReduce(Y, Y);
    nI = I.extend(temp);

    nY = zz_pXi<N+1>(Y);

    nX = zz_pXi<N + 1>(X) - A.a*nI.mulReduce(Zn, nY + nY - Zn);
    sub(nY, nY, Zn);
    cout << "    Time : " << get_time() - start << endl;
}

// WARNING here N is total number of rounds and K is the K'th round processing
template <auto K, auto N>
inline void rec_build_system(zz_pXi<N+1>& nX, zz_pXi<N+1>& nY, ideal<N>& nI, zz_pXi<K>& X, zz_pXi<K>& Y, ideal<K - 1>& I, const anemoi& A)
{
    if constexpr (K == N + 1)
    {
        nX = X;
        nY = Y;
        nI = I;
    }
    else
    {
        cout << "\n[+] Round " << N << endl;
        cout << "[+] degrees : ";
        for (auto& deg : degrees(Y))
            cout << deg << " ";
        cout << endl;
        cout << "[+] Biggest X monom deg : " << deg(getMonomX(Y)) << endl;

        zz_pXi<K + 1> nX2, nY2;
        ideal<K> nI2;
        R_tri(nX2, nY2, nI2, X, Y, I, A);
        rec_build_system<K + 1, N>(nX, nY, nI, nX2, nY2, nI2, A);
    }
}

template <auto N> // N rounds
inline void build_system(zz_pXi<N+1>& nX, zz_pXi<N+1>& nY, ideal<N>& nI)
{
    ideal<0> I0{};
    anemoi A(N);
    zz_pXi<1> X{}; X.Xn();
    zz_pXi<1> Y{}; Y.zero();

    rec_build_system<1, N>(nX, nY, nI, X, Y, I0, A);

    add(nY, nY, nX);
    add(nX, nX, nY);
}


template <auto N>
zz_pX getMonomX(const zz_pXi<N>& Y)
{ return getMonomX(Y.rep[0]); }

template <>
zz_pX getMonomX(const zz_pXi<1>& Y)
{ return Y; }

template <auto N>
inline void reduce_system(zz_pX& nY, zz_pXi<N + 1>& Y, ideal<N>& I)
{
    cout << "\n[+] Round " << N << endl;
    cout << "[+] degrees : ";
    for (auto& deg : degrees(Y))
        cout << deg << " ";
    cout << endl;
    cout << "[+] Biggest X monom deg : " << deg(getMonomX(Y)) << endl;

    double start = get_time();

    if constexpr (N == 1)
    {
        zz_pXi<1> nY2{};
        resultant_P_ideal(nY2, Y, I);
        nY.rep = move(nY2.rep);
    }
    else
    {
        zz_pXi<N> nY2{};
        resultant_P_ideal(nY2, Y, I);
        cout << "    Time : " << get_time() - start << endl;

        reduce_system(nY, nY2, I.base);
    }

}

long heuristic_count = 0;

int main(int argc, char* argv[]){ 
    zz_p::init(P);
    
    anemoi Anemoi(R);
    zz_pXi<R + 1> X, Y;
    ideal<R> I{};

    cout << "[+] Building system" << endl;;
    build_system(X, Y, I);
    cout << "[+] Building system done" << endl << endl;;
    zz_pX eq{};

    cout << "[+] Reducing system" << endl;;
    reduce_system(eq, Y, I);
    cout << "[+] Reducing system done" << endl << endl;;

    cout << "[+] Degree of final poly =" << deg(eq) << endl;
    vec_zz_p roots = zz_pX_roots_opti(eq);

    if (roots.length() == 0) 
    {
        cout << "No solution found for ... (no roots)" << endl;
        return 0;
    }

    zz_p Xa = roots[0];
    zz_p Ya(0);

    cout << endl << "[!] SOLUTION FOUND !!!" << endl;
    cout << "X = " << Xa << endl;
    cout << "Y = " << Ya << endl;

    zz_p Xs(Xa);
    zz_p Ys(Ya);

    Anemoi.encrypt(Xs, Ys);
    cout << "anemoi(" << Xa << ", " << Ya << ", " << R << " rounds) = " << "(" << Xs << ", " << Ys << ")" << endl;

    return 0;
}