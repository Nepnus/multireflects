#include <stdint.h>
#include <math.h>
#include "CubicSplineInterp.h"

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

EXPORT void gen_Teff_distribution(const float* points, uint32_t points_len, 
    float T1, float gdc1, 
    float a1, float b1, float c1,
    float g_abratio1, float g_b2bratio1, float g_cbratio1,
    float* Teff_distribution1){
    
    const float tol = 1.0e-8;

    float T1_m4 = T1 * T1 * T1 * T1;

    uint32_t i;
    float tmp1, tmp2, tmp3;
    float points1_x, points1_y;
    float g_ratio_tmp;

    float g1_a0 = 1.0 / (g_b2bratio1 - 1.0 + tol);
    tmp1 = c1*c1 / (g_abratio1 - 1.0 + tol);
    tmp2 = a1*a1 / (g_cbratio1 - 1.0 + tol);
    tmp3 = a1*a1*c1*c1;
    float g1_a2 = (tmp1 - tmp2) / (tmp3*(a1 + c1)) + g1_a0 * (a1 - c1) / tmp3;
    float g1_a1 = (tmp1 - g1_a0*c1*c1 - g1_a2*a1*tmp3) / tmp3;

    for(i=0; i<points_len; i++){
        points1_x = points[3*i]; points1_y = points[3*i+1];

        tmp1 = points1_y / b1;
        tmp1 *= tmp1;
        tmp2 = points1_x * points1_x;
        g_ratio_tmp = 1.0 + (1.0 - tmp1) / (g1_a0 + g1_a1*tmp2 + g1_a2*tmp2*points1_x);
        tmp1 = T1_m4 * powf(g_ratio_tmp, gdc1);

        Teff_distribution1[i] = logf(tmp1);
    }
}

EXPORT void gen_lc(const float* Teff_distribution1, const float* limb_factor1, const float* areas1, const float* cos_gamma_1, 
    const float* RV1, float bfac1,
    const float* band_a1, const float* band_b1, const float* band_c1, const float* band_d1, float band_hlogT, float band_logT0, 
    const uint32_t* i_ps1, const uint32_t* i_sum1, uint32_t phi_len, float* lc, float* I_distribution_tmp, uint32_t points_len){

    const float light_speed = 299792.458;

    uint32_t i,j,k;
    uint32_t i_sum_single, j_count=0;
    float doppler_factor, I_single, limb_darken, cosgamma_single;
    float tmp1;

    for(i=0; i<points_len; i++)
        I_distribution_tmp[i] = interp_single(Teff_distribution1[i], band_hlogT, band_logT0, band_a1, band_b1, band_c1, band_d1);

    for(i=0; i<phi_len; i++){
        i_sum_single = i_sum1[i];
        doppler_factor = 1.0 - bfac1 * RV1[i] / light_speed;
        lc[i] = 0.0;
        for(j=j_count; j<j_count+i_sum_single; j++){
            k = i_ps1[j];
            cosgamma_single = cos_gamma_1[j];
            tmp1 = sqrt(cosgamma_single);
            limb_darken = 1.0 - limb_factor1[0]*(1.0-tmp1) - limb_factor1[1]*(1.0-cosgamma_single) 
                - limb_factor1[2]*(1.0-tmp1*cosgamma_single) - limb_factor1[3]*(1.0-cosgamma_single*cosgamma_single);
            tmp1 = limb_darken * areas1[k] * cosgamma_single;
            lc[i] += I_distribution_tmp[k] * tmp1;
        }
        lc[i] *= doppler_factor;
        j_count += i_sum_single;
    }
}

