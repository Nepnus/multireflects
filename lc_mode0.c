#include <stdint.h>
#include <math.h>
#include "Fin_interp.h"
#include "CubicSplineInterp.h"

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

float absf_(float x){
    float y = x;
    if(x < 0.0)
        y = -x;
    return y;
}

float minf_(float x, float y){
    float z = x;
    if(y < x)
        z = y;
    return z;
}

EXPORT void gen_cosgamma_multid(float incl, const float* nv, const float* d, uint32_t phi_len,
    float a1, float b1, float b12, float c1,
    float a2, float b2, float b22, float c2,
    const float* unitNormals1, const float* unitNormals2, 
    const float* points1, const float* points2, 
    uint32_t points1_len, uint32_t points2_len,
    float* cos_gamma1, uint32_t* i_ps1, uint32_t* i_sum1,
    float* cos_gamma2, uint32_t* i_ps2, uint32_t* i_sum2){
    
    const float PI = 3.141592653589793;
    const float PI_m2 = PI * 2.0;
    const float PI_d2 = PI / 2.0;
    const float tol = 1.0e-6;
    float incl_rad = incl * PI / 180.0;
    float cos_incl_rad = cosf(incl_rad);
    float sin_incl_rad = sinf(incl_rad);
    float a1_p2 = a1 * a1; float b1_p2 = b1 * b1; float b12_p2 = b12 * b12;
    float a2_p2 = a2 * a2; float b2_p2 = b2 * b2; float b22_p2 = b22 * b22;
    float a1_plus_a2 = a1 + a2;
    uint32_t i_phi, i_points;
    uint32_t ind1=0, ind2=0, ind3, ind4, ind5, count;
    uint8_t eclipse_flag = 0;
    float nv_s, d_s;
    float sin_nv_s, cos_nv_s;
    float nvisual_x, nvisual_y=cos_incl_rad, nvisual_z;
    float pstar_x, pstar_y, pstar_z;
    float tmp1, tmp2, tmp3, tmp4, tmp5;

    for(i_phi=0; i_phi<phi_len; i_phi++){
        nv_s = nv[i_phi]; d_s = d[i_phi];
        sin_nv_s = sinf(nv_s); cos_nv_s = cosf(nv_s);
        nvisual_x = -sin_incl_rad*cos_nv_s;
        nvisual_z = -sin_incl_rad*sin_nv_s;

        // eclipse_judge
        tmp3 = sin_nv_s * sin_nv_s;
        tmp4 = cos_nv_s * cos_nv_s;
        tmp5 = cos_incl_rad * cos_incl_rad;
        tmp4 *= tmp5;
        tmp3 += tmp4;
        tmp3 = sqrtf(tmp3);
        tmp3 *= d_s;
        if(a1_plus_a2 > tmp3){
            if(eclipse_flag == 0){
                tmp3 = minf_(nv_s, PI_m2 - nv_s);
                tmp4 = absf_(nv_s - PI_d2);
                if(tmp3 < tmp4)
                    eclipse_flag = 2;
                else
                    eclipse_flag = 1;
            }
        }
        else{
            eclipse_flag = 0;
        }

        // star1_cosgamma
        count = 0;
        for(i_points=0; i_points<points1_len; i_points++){
            ind3 = 3*i_points; ind4 = 3*i_points+1; ind5 = 3*i_points+2;

            tmp1 = unitNormals1[ind3] * nvisual_x;
            tmp2 = unitNormals1[ind4] * nvisual_y;
            tmp3 = unitNormals1[ind5] * nvisual_z;
            tmp1 += tmp2;
            tmp1 += tmp3;
            if(tmp1 <= 0)
                continue;
            
            if(eclipse_flag == 1){
                pstar_x = points1[ind3]; pstar_y = points1[ind4]; pstar_z = points1[ind5];
                tmp2 = pstar_x - d_s;
                tmp3 = 2.0 * (nvisual_x*tmp2/a2_p2 + nvisual_y*pstar_y/b2_p2 + nvisual_z*pstar_z/b22_p2);
                tmp4 = 4.0 * (nvisual_x*nvisual_x/a2_p2 + nvisual_y*nvisual_y/b2_p2 + nvisual_z*nvisual_z/b22_p2);
                tmp5 = tmp2*tmp2/a2_p2 + pstar_y*pstar_y/b2_p2 + pstar_z*pstar_z/b22_p2 - 1.0;
                tmp2 = tmp3 * tmp3 - tmp4 * tmp5;
                if(tmp2 >= 0.0)
                    continue;
            }

            cos_gamma1[ind1] = tmp1;
            i_ps1[ind1] = i_points;
            count++;
            ind1++;
        }
        i_sum1[i_phi] = count;

        // star2_cosgamma
        count = 0;
        for(i_points=0; i_points<points2_len; i_points++){
            ind3 = 3*i_points; ind4 = 3*i_points+1; ind5 = 3*i_points+2;

            tmp1 = -unitNormals2[ind3] * nvisual_x;
            tmp2 = unitNormals2[ind4] * nvisual_y;
            tmp3 = unitNormals2[ind5] * nvisual_z;
            tmp1 += tmp2;
            tmp1 += tmp3;
            if(tmp1 <= 0)
                continue;
            
            if(eclipse_flag == 2){
                pstar_x = -points2[ind3]; pstar_y = points2[ind4]; pstar_z = points2[ind5];
                tmp2 = pstar_x + d_s;
                tmp3 = 2.0 * (nvisual_x*tmp2/a1_p2 + nvisual_y*pstar_y/b1_p2 + nvisual_z*pstar_z/b12_p2);
                tmp4 = 4.0 * (nvisual_x*nvisual_x/a1_p2 + nvisual_y*nvisual_y/b1_p2 + nvisual_z*nvisual_z/b12_p2);
                tmp5 = tmp2*tmp2/a1_p2 + pstar_y*pstar_y/b1_p2 + pstar_z*pstar_z/b12_p2 - 1.0;
                tmp2 = tmp3 * tmp3 - tmp4 * tmp5;
                if(tmp2 >= 0.0)
                    continue;
            }

            cos_gamma2[ind2] = tmp1;
            i_ps2[ind2] = i_points;
            count++;
            ind2++;
        }
        i_sum2[i_phi] = count;
    }
}

void gen_cosgamma_singled(float incl, const float* nv, float sma, uint32_t phi_len,
    float a1, float b1, float b12, float c1,
    float a2, float b2, float b22, float c2,
    const float* unitNormals1, const float* unitNormals2, 
    const float* points1, const float* points2, 
    uint32_t points1_len, uint32_t points2_len,
    float* cos_gamma1, uint32_t* i_ps1, uint32_t* i_sum1,
    float* cos_gamma2, uint32_t* i_ps2, uint32_t* i_sum2){
    
    const float PI = 3.141592653589793;
    const float PI_m2 = PI * 2.0;
    const float PI_d2 = PI / 2.0;
    const float tol = 1.0e-6;
    float incl_rad = incl * PI / 180.0;
    float cos_incl_rad = cosf(incl_rad);
    float sin_incl_rad = sinf(incl_rad);
    float a1_p2 = a1 * a1; float b1_p2 = b1 * b1; float b12_p2 = b12 * b12;
    float a2_p2 = a2 * a2; float b2_p2 = b2 * b2; float b22_p2 = b22 * b22;
    float a1_plus_a2 = a1 + a2;
    uint32_t i_phi, i_points;
    uint32_t ind1=0, ind2=0, ind3, ind4, ind5, count;
    uint8_t eclipse_flag = 0;
    float nv_s, sin_nv_s, cos_nv_s;
    float nvisual_x, nvisual_y=cos_incl_rad, nvisual_z;
    float pstar_x, pstar_y, pstar_z;
    float tmp1, tmp2, tmp3, tmp4, tmp5;

    for(i_phi=0; i_phi<phi_len; i_phi++){
        nv_s = nv[i_phi];
        sin_nv_s = sinf(nv_s); cos_nv_s = cosf(nv_s);
        nvisual_x = -sin_incl_rad*cos_nv_s;
        nvisual_z = -sin_incl_rad*sin_nv_s;

        // eclipse_judge
        tmp3 = sin_nv_s * sin_nv_s;
        tmp4 = cos_nv_s * cos_nv_s;
        tmp5 = cos_incl_rad * cos_incl_rad;
        tmp4 *= tmp5;
        tmp3 += tmp4;
        tmp3 = sqrtf(tmp3);
        tmp3 *= sma;
        if(a1_plus_a2 > tmp3){
            if(eclipse_flag == 0){
                tmp3 = minf_(nv_s, PI_m2 - nv_s);
                tmp4 = absf_(nv_s - PI_d2);
                if(tmp3 < tmp4)
                    eclipse_flag = 2;
                else
                    eclipse_flag = 1;
            }
        }
        else{
            eclipse_flag = 0;
        }

        // star1_cosgamma
        count = 0;
        for(i_points=0; i_points<points1_len; i_points++){
            ind3 = 3*i_points; ind4 = 3*i_points+1; ind5 = 3*i_points+2;

            tmp1 = unitNormals1[ind3] * nvisual_x;
            tmp2 = unitNormals1[ind4] * nvisual_y;
            tmp3 = unitNormals1[ind5] * nvisual_z;
            tmp1 += tmp2;
            tmp1 += tmp3;
            if(tmp1 <= 0)
                continue;
            
            if(eclipse_flag == 1){
                pstar_x = points1[ind3]; pstar_y = points1[ind4]; pstar_z = points1[ind5];
                tmp2 = pstar_x - sma;
                tmp3 = 2.0 * (nvisual_x*tmp2/a2_p2 + nvisual_y*pstar_y/b2_p2 + nvisual_z*pstar_z/b22_p2);
                tmp4 = 4.0 * (nvisual_x*nvisual_x/a2_p2 + nvisual_y*nvisual_y/b2_p2 + nvisual_z*nvisual_z/b22_p2);
                tmp5 = tmp2*tmp2/a2_p2 + pstar_y*pstar_y/b2_p2 + pstar_z*pstar_z/b22_p2 - 1.0;
                tmp2 = tmp3 * tmp3 - tmp4 * tmp5;
                if(tmp2 >= 0.0)
                    continue;
            }

            cos_gamma1[ind1] = tmp1;
            i_ps1[ind1] = i_points;
            count++;
            ind1++;
        }
        i_sum1[i_phi] = count;

        // star2_cosgamma
        count = 0;
        for(i_points=0; i_points<points2_len; i_points++){
            ind3 = 3*i_points; ind4 = 3*i_points+1; ind5 = 3*i_points+2;

            tmp1 = -unitNormals2[ind3] * nvisual_x;
            tmp2 = unitNormals2[ind4] * nvisual_y;
            tmp3 = unitNormals2[ind5] * nvisual_z;
            tmp1 += tmp2;
            tmp1 += tmp3;
            if(tmp1 <= 0)
                continue;
            
            if(eclipse_flag == 2){
                pstar_x = -points2[ind3]; pstar_y = points2[ind4]; pstar_z = points2[ind5];
                tmp2 = pstar_x + sma;
                tmp3 = 2.0 * (nvisual_x*tmp2/a1_p2 + nvisual_y*pstar_y/b1_p2 + nvisual_z*pstar_z/b12_p2);
                tmp4 = 4.0 * (nvisual_x*nvisual_x/a1_p2 + nvisual_y*nvisual_y/b1_p2 + nvisual_z*nvisual_z/b12_p2);
                tmp5 = tmp2*tmp2/a1_p2 + pstar_y*pstar_y/b1_p2 + pstar_z*pstar_z/b12_p2 - 1.0;
                tmp2 = tmp3 * tmp3 - tmp4 * tmp5;
                if(tmp2 >= 0.0)
                    continue;
            }

            cos_gamma2[ind2] = tmp1;
            i_ps2[ind2] = i_points;
            count++;
            ind2++;
        }
        i_sum2[i_phi] = count;
    }
}

EXPORT void gen_Teff_distribution_multid(const float* points, const float* thetas, uint32_t points_len, 
    float T1, float T2, float A1, float gdc1, 
    float a1, float b1, float c1,
    float g_abratio1, float g_b2bratio1, float g_cbratio1,
    const float* d, uint32_t d_len, const uint32_t* i_ps1, const uint32_t* i_sum1,
    float Fin1_d_h, float Fin1_d_0, float Fin1_theta2_h, float Fin1_theta2_max, uint32_t Fin1_d_num, uint32_t Fin1_theta2_num,
    const float* Fin1_A, const float* Fin1_B, const float* Fin1_C, const float* Fin1_D,
    float* Fin1_a, float* Fin1_b, float* Fin1_c, float* Fin1_d,
    float* Teff_distribution1, float* incident_light_factor1, float* Teff_distribution_tmp
    ){
    
    const float tol = 1.0e-8;
    const float PI = 3.141592653589793;

    float factor1 = (1.0 - A1) / PI * T2 * T2 * T2 * T2;
    float T1_m4 = T1 * T1 * T1 * T1;

    uint32_t i,j,k;
    uint32_t i_sum_single, j_count=0;
    float tmp1, tmp2, tmp3;
    float points1_x, points1_y;
    float incident_factor, factor2;
    float g_ratio_tmp;

    float g1_a0 = 1.0 / (g_b2bratio1 - 1.0 + tol);
    tmp1 = c1*c1 / (g_abratio1 - 1.0 + tol);
    tmp2 = a1*a1 / (g_cbratio1 - 1.0 + tol);
    tmp3 = a1*a1*c1*c1;
    float g1_a2 = (tmp1 - tmp2) / (tmp3*(a1 + c1)) + g1_a0 * (a1 - c1) / tmp3;
    float g1_a1 = (tmp1 - g1_a0*c1*c1 - g1_a2*a1*tmp3) / tmp3;

    for(i=0; i<points_len; i++)
        Teff_distribution_tmp[i] = -1.0;

    for(i=0; i<d_len; i++){
        i_sum_single = i_sum1[i];
        ddirection_interp(d[i], Fin1_d_h, Fin1_d_0, Fin1_theta2_num, Fin1_d_num, Fin1_A, Fin1_B, Fin1_C, Fin1_D, Fin1_a);
        get_thetadirection_factors(Fin1_a, Fin1_theta2_h, Fin1_theta2_num, Fin1_b, Fin1_c, Fin1_d);
        for(j=j_count; j<j_count+i_sum_single; j++){
            k = i_ps1[j];
            factor2 = Teff_distribution_tmp[k];

            if(factor2 < 0.0){
                points1_x = points[3*k]; points1_y = points[3*k+1];

                tmp1 = points1_y / b1;
                tmp1 *= tmp1;
                tmp2 = points1_x * points1_x;
                g_ratio_tmp = 1.0 + (1.0 - tmp1) / (g1_a0 + g1_a1*tmp2 + g1_a2*tmp2*points1_x);
                factor2 = T1_m4 * powf(g_ratio_tmp, gdc1);
                Teff_distribution_tmp[k] = factor2;
            }

            incident_factor = thetadirection_interp(thetas[k], Fin1_theta2_h, Fin1_theta2_max, Fin1_a, Fin1_b, Fin1_c, Fin1_d);
            Teff_distribution1[j] = logf(incident_factor * factor1 + factor2);
            incident_light_factor1[j] = incident_factor;
        }
        j_count += i_sum_single;
    }
}

EXPORT void gen_Teff_distribution_singled(const float* points, const float* thetas, uint32_t points_len, 
    float T1, float T2, float A1, float gdc1, 
    float a1, float b1, float c1,
    float g_abratio1, float g_b2bratio1, float g_cbratio1,
    float Fin1_theta2_h, float Fin1_theta2_max,
    const float* Fin1_a, const float* Fin1_b, const float* Fin1_c, const float* Fin1_d,
    float* Teff_distribution1, float* incident_light_factor1){

    const float tol = 1.0e-8;
    const float PI = 3.141592653589793;

    float factor1 = (1.0 - A1) / PI * T2 * T2 * T2 * T2;
    float T1_m4 = T1 * T1 * T1 * T1;

    uint32_t i;
    float tmp1, tmp2, tmp3;
    float points1_x, points1_y;
    float incident_factor, factor2;
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
        factor2 = T1_m4 * powf(g_ratio_tmp, gdc1);

        incident_factor = thetadirection_interp(thetas[i], Fin1_theta2_h, Fin1_theta2_max, Fin1_a, Fin1_b, Fin1_c, Fin1_d);
        Teff_distribution1[i] = logf(incident_factor * factor1 + factor2);
        incident_light_factor1[i] = incident_factor;
    }
}

EXPORT void gen_lc_multid(const float* Teff_distribution1, const float* incident_factor1, const float* limb_factor1, const float* areas1, const float* cos_gamma_1, 
    const float* RV1, float bfac1, float A1, float T0,
    const float* band_a1, const float* band_b1, const float* band_c1, const float* band_d1, float band_hlogT, float band_logT0, 
    const float* band_a0, const float* band_b0, const float* band_c0, const float* band_d0,
    const uint32_t* i_ps1, const uint32_t* i_sum1, uint32_t phi_len, float* lc){

    const float light_speed = 299792.458;
    const float PI = 3.141592653589793;
    float outgoing_light = interp_single(4.0*logf(T0), band_hlogT, band_logT0, band_a0, band_b0, band_c0, band_d0);
    outgoing_light *= A1 / PI / 2.0;

    uint32_t i,j,k;
    uint32_t i_sum_single, j_count=0;
    float doppler_factor, I_single, limb_darken, cosgamma_single;
    float tmp1;

    for(i=0; i<phi_len; i++){
        i_sum_single = i_sum1[i];
        doppler_factor = 1.0 - bfac1 * RV1[i] / light_speed;
        lc[i] = 0.0;
        for(j=j_count; j<j_count+i_sum_single; j++){
            k = i_ps1[j];
            cosgamma_single = cos_gamma_1[j];
            I_single = interp_single(Teff_distribution1[j], band_hlogT, band_logT0, band_a1, band_b1, band_c1, band_d1);
            tmp1 = sqrt(cosgamma_single);
            limb_darken = 1.0 - limb_factor1[0]*(1.0-tmp1) - limb_factor1[1]*(1.0-cosgamma_single) 
                - limb_factor1[2]*(1.0-tmp1*cosgamma_single) - limb_factor1[3]*(1.0-cosgamma_single*cosgamma_single);
            tmp1 = limb_darken * areas1[k] * cosgamma_single;
            lc[i] += I_single * tmp1;
            lc[i] += outgoing_light * incident_factor1[j] * tmp1;
        }
        lc[i] *= doppler_factor;
        j_count += i_sum_single;
    }
}

EXPORT void gen_lc_singled(const float* Teff_distribution1, const float* incident_factor1, const float* limb_factor1, const float* areas1, const float* cos_gamma_1, 
    const float* RV1, float bfac1, float A1, float T0,
    const float* band_a1, const float* band_b1, const float* band_c1, const float* band_d1, float band_hlogT, float band_logT0, 
    const float* band_a0, const float* band_b0, const float* band_c0, const float* band_d0,
    const uint32_t* i_ps1, const uint32_t* i_sum1, uint32_t phi_len, float* lc, float* I_distribution_tmp, uint32_t points_len){

    const float light_speed = 299792.458;
    const float PI = 3.141592653589793;
    float outgoing_light = interp_single(4.0*logf(T0), band_hlogT, band_logT0, band_a0, band_b0, band_c0, band_d0);
    outgoing_light *= A1 / PI / 2.0;

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
            lc[i] += outgoing_light * incident_factor1[k] * tmp1;
        }
        lc[i] *= doppler_factor;
        j_count += i_sum_single;
    }
}






