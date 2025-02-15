#include <iostream>
#include <fstream>
#include <NTL/ALL_FEATURES.h>
#include <NTL/lzz_p.h>
#include <NTL/vec_vec_lzz_p.h>
#include "lzz_pXi_triangular_ideal.h"
#include "lzz_pX_opti_roots.h"
#include "lzz_pXi.h"
#include "special_resultant.h"
#include "griffin.h"
#include "griffin_mat.h"
#include "pp.h"
#include "vec_lzz_pXi.h"
#include "mat_lzz_pXi.h"
#include <unistd.h>

extern Vec<zz_p> BYPASS_MUL_GRIFFIN;
extern Vec<zz_p> BYPASS_ADD_GRIFFIN;

template <auto N>
using ideal = zz_pXi_triangular_ideal<GRIFFIN_D, N>;

NTL_CLIENT

template <auto N>
inline zz_pXi<N+1> reduce_G_i(const zz_pXi<N+1>& x, const zz_pXi<N+1>& y, zz_pXi<N>& z, int i, const ideal<N - 1>& I, const ideal<N>& nI,
    const zz_pXi<N + 1>& T1, const zz_pXi<N + 1>& T2, const zz_pXi<N + 1>& T3)
{
    zz_pXi<N+1> T4 = GAMMAS_GRIFFIN[i]*x + y;
    zz_pXi<N+1> z_up(z);
    zz_pXi<N+1> T5 = nI.mulReduce(T4, z_up);
    
    zz_pXi<N> red_z_sq_low = I.sqrReduce(z);
    zz_pXi<N+1> red_z_sq(red_z_sq_low);
    return GAMMAS_GRIFFIN[i] * GAMMAS_GRIFFIN[i] * T1 + 2 * GAMMAS_GRIFFIN[i] * T3 + ALPHAS_GRIFFIN[i]*(T4 + z_up) + zz_p(2)*T5 + T2 +
        red_z_sq + BETAS_GRIFFIN[i];
}

template <auto N>
inline void R_tri(Vec<zz_pXi<N+1>>& nState, ideal<N>& nI, Vec<zz_pXi<N>>& state, ideal<N - 1>& I, const griffin& G)
{   
    double start;

    // csts
    cout << " - non-linear of round " << N + GRIFFIN_BYPASS << endl;
    cout << "   [+] Non-linear (not the SBOX)" << endl;
    start = get_time();

    cout << "   [+] Non-linear wire 0 (actually the SBOX)" << endl;
    nState[0].Xn();

    cout << "   [+] Non-linear wire 1" << endl;
    zz_pXi<N> tmp1 = I.mulReduce(I.mulReduce(state[1], state[1]), state[1]);
    nState[1] = tmp1.extend();
    

    // notation : 
    zz_pXi<N + 1> y0 = nState[0];
    zz_pXi<N + 1> y1 = nState[1];

    nI = I.extend(state[0]);

    cout << "   [+] Precomputation for the other wires" << endl;
    zz_pXi<N + 1> T1 = nI.sqrReduce(y0);

    // this is done as an effort to not reduce in an "unnecessary big ideal"
    zz_pXi<N> T2_low = I.sqrReduce(tmp1);
    zz_pXi<N + 1> T2 = T2_low.extend();
    zz_pXi<N + 1> T3 = nI.mulReduce(y0, y1);

    // now expression should contain new variables
    zz_pXi<N> zero{}; zero.zero();
    cout << "   [+] Non-linear wire 2" << endl;
    nState[2] = nI.mulReduce(state[2].extend(), reduce_G_i(y0, y1, zero, 2, I, nI, T1, T2, T3));

    for (int i = 3; i < GRIFFIN_B; i++)
    {
        cout << "   [+] Non-linear wire " << i << endl;
        nState[i] = nI.mulReduce(state[i].extend(), reduce_G_i(y0, y1, state[i - 1], i, I, nI, T1, T2, T3));
    }

    cout << "    Time : " << get_time() - start << endl;

    // linear
    cout << " - Linear of round " << N + GRIFFIN_BYPASS << endl;
    start = get_time();
    G.linear_layer(nState, N - 1 + GRIFFIN_BYPASS);
    cout << "    Time : " << get_time() - start << endl;
}

// WARNING here N is total number of rounds and K is the K'th round processing
template <auto K, auto N>
inline void rec_build_system(Vec<zz_pXi<N+1>>& nState, ideal<N>& nI, Vec<zz_pXi<K>>& state, ideal<K - 1>& I, const griffin& G)
{
    if constexpr (K == N + 1)
    {
        nState = state;
        nI = I;
    }
    else
    {
        cout << "\n[+] Round " << K + GRIFFIN_BYPASS << endl;

        Vec<zz_pXi<K + 1>> nState2{};
        nState2.SetLength(GRIFFIN_B);

        ideal<K> nI2;
        R_tri(nState2, nI2, state, I, G);
        rec_build_system<K + 1, N>(nState, nI, nState2, nI2, G);
    }
}

template <auto N> // N rounds
inline void build_system(Vec<zz_pXi<N+1>>& output, ideal<N>& nI)
{
    ideal<0> I0{};
    griffin G{};

    zz_pXi<1> X{}; X.Xn();

    Vec<zz_pXi<1>> state;
    state.SetLength(GRIFFIN_B);
    for (int i = 0; i < GRIFFIN_B; i++)
        state[i] = BYPASS_MUL_GRIFFIN[i] * X + BYPASS_ADD_GRIFFIN[i];

    std::cout << "input : " << pp(state) << std::endl << std::endl;
 
    // the bypassable part
    mul(state, G.MDS, state);
    for (int i = 0; i < GRIFFIN_BYPASS; ++i)
    {
        cout << " - SBOX of round " << i + 1 << endl;
        G.s_box(state);
        cout << " - Linear of round " << i + 1 << endl;
        G.linear_layer(state, i);
    }
    rec_build_system<1, N>(output, nI, state, I0, G);
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
    cout << "\n[+] Round " << N + GRIFFIN_BYPASS << endl;
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

int main(int argc, char* argv[])
{
    zz_p::init(GRIFFIN_P);
    init_csts_griffin();

    // attack
    Vec<zz_pXi<GRIFFIN_R + 1 - GRIFFIN_BYPASS>> state{};
    state.SetLength(GRIFFIN_B);

    ideal<GRIFFIN_R - GRIFFIN_BYPASS> I{};

    cout << "[+] Building system" << endl;;
    build_system(state, I);
    cout << "[+] Building system done" << endl << endl;;
    zz_pX eq{};

    cout << "[+] Reducing system" << endl;;
    reduce_system(eq, state[0], I);
    cout << "[+] Reducing system done" << endl << endl;;

    cout << "[+] degree of univariate polynomial is : " << deg(eq) << endl;

    vec_zz_p roots = zz_pX_roots_opti(eq);

    if (roots.length() == 0) 
    {
        cout << "No solution found for ... (no roots)" << endl;
        return 0;
    }

    for (auto &&r : roots)
    {
        cout << endl << "[!] SOLUTION FOUND !!!" << endl;
        cout << "X = " << r << endl;

        Vec<zz_p> stateVerif{};
        stateVerif.SetLength(GRIFFIN_B);
        for (int i = 0; i < GRIFFIN_B; i++)
            stateVerif[i] = BYPASS_MUL_GRIFFIN[i] * r + BYPASS_ADD_GRIFFIN[i];

        cout << "INPUT :" << stateVerif << endl;
        griffin G{};
        G.encrypt(stateVerif);
        cout << "OUTPUT :" << stateVerif << endl << endl;
    }
}