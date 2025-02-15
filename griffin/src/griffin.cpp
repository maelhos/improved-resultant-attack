#include "griffin.h"
#include <assert.h>

Vec<Vec<zz_p>> CSTS_GRIFFIN;
Vec<zz_p> ALPHAS_GRIFFIN;
Vec<zz_p> BETAS_GRIFFIN;
Vec<zz_p> GAMMAS_GRIFFIN;
Vec<zz_p> BYPASS_MUL_GRIFFIN;
Vec<zz_p> BYPASS_ADD_GRIFFIN;

griffin::griffin() : MDS(INIT_SIZE, GRIFFIN_B, GRIFFIN_B) 
{
    zz_pPush push(zz_p::modulus() - 1);
    inv_d = (zz_p(1) / zz_p(GRIFFIN_D))._zz_p__rep;
    build_mds();
}

void griffin::build_mds()
{
    if constexpr (GRIFFIN_B == 3)
    {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                MDS[i][j] = i == j ? 2 : 1;
    }
    else
    {   
        // first build the matrix M'
        Mat<zz_p> Mp(INIT_SIZE, GRIFFIN_B, GRIFFIN_B);
        int tp = GRIFFIN_B / 4;

        int raw_M4[16] = {5, 7, 1, 3, 4, 6, 1, 1, 1, 3, 5, 7, 1, 1, 4, 6};
        for (int b = 0; b < tp; ++b)
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    Mp[4*b + i][4*b + j] = raw_M4[4*i+j];

        if constexpr (GRIFFIN_B == 4)
            MDS = Mp;
        else // then t = 4*t' >= 8
        {
            Mat<zz_p> Mpp(INIT_SIZE, GRIFFIN_B, GRIFFIN_B);

            // build circ([2*I, I, ..., I]), Mael : I did not make this
            for (int b1 = 0; b1 < tp; ++b1)
                for (int b2 = 0; b2 < tp; ++b2)
                    for (int i = 0; i < 4; ++i)
                        for (int j = 0; j < 4; ++j)
                            Mpp[4*b1+i][4*b2+j] = b1 == b2 ? (i == j ? 2 : 0) : (i == j ? 1 : 0);
            
            MDS = Mp * Mpp;
        }   
    }
}

void griffin::encrypt(griffin_state& s){ 
    Vec<zz_p> res(s);
    mul(s, MDS, s);
    for (int i = 0; i < GRIFFIN_R; ++i)
    {
        s_box(s);
        linear_layer(s, i);
    }
}