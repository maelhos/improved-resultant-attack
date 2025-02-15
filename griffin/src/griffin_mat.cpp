#include "griffin_mat.h"

void init_csts_griffin(){
    stringstream ss;

    ss << GRIFFIN_CSTS;
    ss >> CSTS_GRIFFIN;

    ss << GRIFFIN_ALPHAS;
    ss >> ALPHAS_GRIFFIN;

    ss << GRIFFIN_BETAS;
    ss >> BETAS_GRIFFIN;

    ss << GRIFFIN_GAMMAS;
    ss >> GAMMAS_GRIFFIN;

    ss << GRIFFIN_BYPASS_MUL;
    ss >> BYPASS_MUL_GRIFFIN;

    ss << GRIFFIN_BYPASS_ADD;
    ss >> BYPASS_ADD_GRIFFIN;
}