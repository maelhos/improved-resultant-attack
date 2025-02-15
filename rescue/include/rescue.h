
#ifndef RESCUE__H
#define RESCUE__H

#include <NTL/lzz_p.h>
#include <NTL/mat_lzz_p.h>
#include <array>
#include "rescue_mat.h"

NTL_CLIENT


using rescue_state = Vec<zz_p>;

class rescue
{
    long inv_d;
public:
    long round_constants_counter;

    rescue();
    void encrypt(rescue_state& s);
    void round(rescue_state& s);
    ~rescue() {}
};

#endif