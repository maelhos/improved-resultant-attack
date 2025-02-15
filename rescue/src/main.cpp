
#include "rescue.h"
#include "rescue_mat.h"
#include <unistd.h>
#include <thread>
#include "lzz_pXi_triangular_ideal.h"
#include "mat_lzz_pXi.h"
#include "lzz_pX_opti_roots.h"
#include "special_resultant.h"
#include "lzz_pXi.h"
#include "pp.h"
#include "vec_lzz_pXi.h"
#include "mat_lzz_pXi.h"

template <auto N>
using ideal = zz_pXi_triangular_ideal<RESCUE_D, N>;

long ctr_cst;

template <auto N>
inline void R_tri(Vec<zz_pXi<N+3>>& nState, ideal<N+2>& nI, Vec<zz_pXi<N>>& state, ideal<N - 1>& I)
{   
    Vec<zz_pXi<N>> nstate_p{};
    nstate_p.SetLength(RESCUE_T);

    cout << " - SBOX of round " << N << endl;
    for (long i = 0; i < RESCUE_T; i++)
    {
        I.sqrReduce(nstate_p[i], state[i]);
        I.mulReduce(nstate_p[i], nstate_p[i], state[i]);
    }

    cout << " - Fisrt Linear of round " << N << endl;
    mul(nstate_p, MDS, nstate_p);

    cout << " - First Constant of round " << N << endl;
    for (long i = 0; i < RESCUE_T; i++)
        nstate_p[i] += round_constants[ctr_cst++];

    // second half

    cout << " - INV-SBOX of round "  << N << endl;
    zz_pXi<N+1> tmp1{}; tmp1.Xn();
    ideal<N> tmpId1 = I.extend(nstate_p[0]);

    zz_pXi<N+2> tmp2{}; tmp2.Xn();
    ideal<N+1> tmpId2 = tmpId1.extend(nstate_p[1].extend());

    nState[0] = tmp1.extend().extend();
    nState[1] = tmp2.extend();
    nState[2].Xn();
    nI = tmpId2.extend(nstate_p[2].extend().extend());

    cout << " - Second Linear of round " << N << endl;
    mul(nState, MDS, nState);

    cout << " - Second Constant of round " << N << endl;
    for (long i = 0; i < RESCUE_T; i++)
        nState[i] += round_constants[ctr_cst++];
}

// WARNING here N is total number of rounds and K is the K'th round processing
template <auto K, auto N>
inline void rec_build_system(Vec<zz_pXi<N+1>>& nState, ideal<N>& nI, Vec<zz_pXi<K>>& state, ideal<K - 1>& I)
{
    if constexpr (K == N + 1)
    {
        nState = state;
        nI = I;
    }
    else
    {
        cout << "\n[+] Round " << K << endl;

        Vec<zz_pXi<K + 3>> nState2{};
        nState2.SetLength(RESCUE_T);

        ideal<K + 2> nI2;
        R_tri(nState2, nI2, state, I);
        rec_build_system<K + 3, N>(nState, nI, nState2, nI2);
    }
}

template <auto N> // N rounds
inline void build_system(Vec<zz_pXi<N+1>>& output, ideal<N>& nI)
{
    ideal<0> I0{};
    ctr_cst = 2 * RESCUE_T; // bypass 1 one round
    zz_pXi<1> X{}; X.Xn();

    Vec<zz_pXi<1>> state;
    state.SetLength(RESCUE_T);

    // for now simple
    state[0] = RESCUE_ALPHA0 * X + RESCUE_BETA0;
    state[1] = RESCUE_ALPHA1 * X + RESCUE_BETA1;
    state[2] = RESCUE_ALPHA2 * X + RESCUE_BETA2;

    std::cout << "input : " << pp(state) << std::endl << std::endl;
    rec_build_system<1, N>(output, nI, state, I0);
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
    cout << "\n[+] Round " << N  << endl;
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

int main(int argc, char** argv){
    zz_p::init(RESCUE_P);
    init_csts_rescue();

    // attack
    Vec<zz_pXi<3*(RESCUE_R) + 1>> state{};
    state.SetLength(RESCUE_T);

    ideal<3*RESCUE_R> I{};

    cout << "[+] Building system" << endl;;
    build_system(state, I);
    cout << "[+] Building system done" << endl << endl;;
    zz_pX eq{};

    cout << "[+] Reducing system" << endl;
    cout << pp(state[2]) << endl;
    reduce_system(eq, state[2], I);
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
    }

    cout << "[+] Heuristic count : " << heuristic_count << endl;
    return 0;
}