#include "CubicSplineInterp.h"
#include "Fin_interp.h"
#include <stdint.h>

void get_ddirection_factors(const float* factors, uint32_t theta2_num, uint32_t d_num, float d_h, float* A, float* B, float* C, float* D){
    uint32_t theta_i, A_ind;

    for(theta_i=0; theta_i<theta2_num; theta_i++){
        A_ind = theta_i * (d_num - 1);
        getInterpFactors_bycol(&factors[theta_i], d_h, d_num, theta2_num, &A[A_ind], &B[A_ind], &C[A_ind], &D[A_ind]);
    }

}

void ddirection_interp(float d_interp, float d_h, float d_0, uint32_t theta2_num, uint32_t d_num,
    const float* A, const float* B, const float* C, const float* D, float* y){
    uint32_t theta_i, A_ind;

    for(theta_i=0; theta_i<theta2_num; theta_i++){
        A_ind = theta_i * (d_num - 1);
        y[theta_i] = interp_single(d_interp, d_h, d_0, &A[A_ind], &B[A_ind], &C[A_ind], &D[A_ind]);
    }

}

void get_thetadirection_factors(const float* factors, float theta_h, uint32_t theta2_num, float* b, float* c, float* d){
    getInterpFactors(factors, theta_h, theta2_num, b, c, d);
}

float thetadirection_interp(float theta_interp, float theta_h, float theta_max,
    const float* a, const float* b, const float* c, const float* d){
    float y = 0.0;
    if(theta_interp < theta_max)
        y = interp_single(theta_interp, theta_h, 0.0, a, b, c, d);
    return y;
}



