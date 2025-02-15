
#ifndef GRIFFIN__H
#define GRIFFIN__H

#include <NTL/lzz_p.h>
#include <NTL/mat_lzz_p.h>
#include <array>
#include "griffin_mat.h"
#include "lzz_pXi.h"
#include "vec_lzz_pXi.h"
#include "mat_lzz_pXi.h"
#include "pp.h"

NTL_CLIENT


using griffin_state = Vec<zz_p>;

class griffin
{
public:
    Mat<zz_p> MDS;
    long inv_d;
    void build_mds();

    griffin();
    void encrypt(griffin_state& s);

    template <auto N, typename Y>
    zz_pXi<N> G_i(const zz_pXi<N>& x, const zz_pXi<N>& y, const Y& z, int i) const
    {
        zz_pXi<N> ret;
        mul<N>(ret, GAMMAS_GRIFFIN[i]*x + y + z, 
            GAMMAS_GRIFFIN[i]*x + y + z + ALPHAS_GRIFFIN[i]);
        return ret + BETAS_GRIFFIN[i];
    }

    template <typename Y>
    zz_p G_i(const zz_p& x, const zz_p y, const Y& z, int i) const
    {
        return (GAMMAS_GRIFFIN[i]*x + y + z)*(GAMMAS_GRIFFIN[i]*x + y + z + ALPHAS_GRIFFIN[i]) + BETAS_GRIFFIN[i];
    }

    template <typename T>
    void s_box(Vec<T>& s) const
    {
        Vec<T> x(s);
        s[0] = power(s[0], inv_d);
        s[1] = power(s[1], GRIFFIN_D);
        s[2] *= G_i(s[0], s[1], 0, 2);
        for (int i = 3; i < GRIFFIN_B; i++)
            s[i] *= G_i(s[0], s[1], x[i - 1], i);
    }

    template <typename T>
    void linear_layer(Vec<T>& s, long i) const
    {
        mul(s, MDS, s);
        if (i != GRIFFIN_R - 1) 
            add(s, s, CSTS_GRIFFIN[i]);
    }

    ~griffin() {}
};

#endif