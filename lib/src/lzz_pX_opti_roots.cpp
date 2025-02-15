#include "lzz_pX_opti_roots.h"

void PowerXMod2(zz_pX& hh, const long& e, const zz_pXModulus& F)
{
    if (F.n < 0) LogicError("PowerXMod: uninitialized modulus");
    cout << "fft = " << F.UseFFT << endl;
    if (e == 0) {
        set(hh);
        return;
    }

    long n = NumBits(e);
    long i;

    zz_pX h, h1;

    h.SetMaxLength(F.n);
    set(h);

    auto total_steps = n;
    size_t steps_completed = 0;
    double startime = get_time();

    for (i = n - 1; i >= 0; i--) {
        if (bit(e, i)) {
            SqrMod(h1, h, F);
            MulByXMod(h, h1, F);
        }
        else
            SqrMod(h, h, F);

        if (deg(h) < F.n - 3)
        {
            cout << "Still small polys... " << get_time() - startime << "s" << endl;
            startime = get_time();
            total_steps--;
            continue;
        }
        
        steps_completed++;
        double elapsed = get_time() - startime;
        double time_estim_total = (total_steps * elapsed) / steps_completed;
        double time_estim_rest = time_estim_total - elapsed;

        std::cout << "\rProgress: " << steps_completed << " of " << total_steps << " (" << std::fixed << std::setprecision(1) << (100.0*steps_completed/total_steps) << "%) time elapsed : " << elapsed << "s, " << time_estim_rest << "s left / total " << time_estim_total << "s" << std::flush;

    }
    cout << endl << flush;
    
    if (e < 0) InvMod(h, h, F);

    hh = h;
}

vec_zz_p zz_pX_roots_opti(const zz_pX& P){
    zz_pX Q_;
    SetCoeff(Q_, 1, to_zz_p(1));

    zz_pX Q;
    double start = get_time();
    PowerXMod2(Q, zz_p::modulus(), P);
    cout << "Powmod Time : " << get_time() - start << endl;
    Q -= Q_;

    zz_pX RR;
    start = get_time();
    GCD(RR, P, Q);
    cout << "GCD Time : " << get_time() - start << endl;

    vec_zz_p roots;
    start = get_time();
    FindRoots(roots, RR);
    cout << "roots Time : " << get_time() - start << endl;

    return roots;
}