#include "rescue_mat.h"

void init_csts_rescue(){
    stringstream ss;

    ss << MDSS;
    ss >> MDS;

    ss << CSTS;
    ss >> round_constants;
}