#include <stdint.h>
#include <math.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

EXPORT void tensorcore(const float* points1, const float* areas1, const float* unitNormals1, const float* thetas1,
    const float* P2, const float* np2, const float* d, uint32_t p1len, uint32_t theta2len, uint32_t dlen, 
    float max_theta1, float u1, float gdc1, float a1, float b1, float c1, float g_ratio1, float g_ratio2, float g_ratio3, float* factor){

    const float tol = 1.0e-8;

    uint32_t i,j,k;
    uint32_t i_ind0, i_ind1, i_ind2, j_ind0, j_ind1;
    uint8_t flag;
    float tmp1, tmp2, tmp3;
    float points1_x, points1_y, points1_z;
    float unitNormals1_x, unitNormals1_y, unitNormals1_z;
    float areas1_i;
    float P2_x, P2_y, P2_z;
    float np2_x, np2_y, np2_z;
    float vector_P1P2_x, vector_P1P2_y, vector_P1P2_z, len_P1P2;
    float cos_alpha, cos_beta;
    float limbdarken, factor_single, gravitydarken;

    float g1_a0 = 1.0 / (g_ratio1 - 1.0 + tol);
    tmp1 = c1*c1 / (g_ratio2 - 1.0 + tol);
    tmp2 = a1*a1 / (g_ratio3 - 1.0 + tol);
    tmp3 = a1*a1*c1*c1;
    float g1_a2 = (tmp1 - tmp2) / (tmp3*(a1 + c1)) + g1_a0 * (a1 - c1) / tmp3;
    float g1_a1 = (tmp1 - g1_a0*c1*c1 - g1_a2*a1*tmp3) / tmp3;

    for(i=0; i<p1len; i++){
        if(thetas1[i] >= max_theta1)
            continue;
        i_ind2 = i*3+2;
        points1_z = points1[i_ind2];
        if(points1_z < 0.0)
            continue;
        i_ind0 = i*3; i_ind1 = i_ind0+1;
        points1_x = points1[i_ind0]; points1_y = points1[i_ind1]; 
        unitNormals1_x = unitNormals1[i_ind0]; unitNormals1_y = unitNormals1[i_ind1]; unitNormals1_z = unitNormals1[i_ind2];
        areas1_i = areas1[i];

        flag = 1;

        for(j=0; j<theta2len; j++){
            j_ind0 = j*2; j_ind1 = j_ind0 + 1;
            P2_x = P2[j_ind0]; P2_y = P2[j_ind1];
            np2_x = np2[j_ind0]; np2_y = np2[j_ind1];
            for(k=0; k<dlen; k++){
                vector_P1P2_x = P2_x - points1_x + d[k];
                vector_P1P2_y = P2_y - points1_y;
                vector_P1P2_z = - points1_z;

                tmp1 = vector_P1P2_x * unitNormals1_x;
                tmp2 = vector_P1P2_y * unitNormals1_y;
                tmp3 = vector_P1P2_z * unitNormals1_z;
                tmp1 += tmp2;
                cos_alpha = tmp1 + tmp3;
                if(cos_alpha <= 0.0)
                    continue;
                
                tmp1 = vector_P1P2_x * np2_x;
                tmp2 = vector_P1P2_y * np2_y;
                cos_beta = tmp1 + tmp2;
                if(cos_beta >= 0.0)
                    continue;

                if(flag){
                    flag = 0;

                    tmp1 = points1_y / b1;
                    tmp1 *= tmp1;
                    tmp2 = points1_x * points1_x;
                    gravitydarken = 1.0 + (1.0 - tmp1) / (g1_a0 + g1_a1*tmp2 + g1_a2*tmp2*points1_x);
                    gravitydarken = powf(gravitydarken, gdc1);
                }

                tmp1 = vector_P1P2_x * vector_P1P2_x;
                tmp2 = vector_P1P2_y * vector_P1P2_y;
                tmp3 = vector_P1P2_z * vector_P1P2_z;
                tmp1 += tmp2;
                tmp1 += tmp3;
                len_P1P2 = sqrtf(tmp1);
                cos_alpha /= len_P1P2;
                cos_beta /= len_P1P2;
                cos_beta *= -1.0;

                limbdarken = 1.0 - u1 * (1.0 - cos_alpha);
                factor_single = limbdarken * gravitydarken;
                factor_single *= cos_alpha;
                factor_single *= cos_beta;
                factor_single *= areas1_i;
                factor_single /= len_P1P2;
                factor_single /= len_P1P2;
                factor[k*theta2len+j] += factor_single;
            }
        }
    }
}

EXPORT void matrixcore(const float* points1, const float* areas1, const float* unitNormals1, const float* thetas1,
    const float* P2, const float* np2, uint32_t p1len, uint32_t theta2len,
    float max_theta1, float sma, float u1, float gdc1, float a1, float b1, float c1, float g_ratio1, float g_ratio2, float g_ratio3, float* factor){

    const float tol = 1.0e-8;

    uint32_t i,j;
    uint32_t i_ind0, i_ind1, i_ind2, j_ind0, j_ind1;
    uint8_t flag;
    float tmp1, tmp2, tmp3;
    float points1_x, points1_y, points1_z;
    float unitNormals1_x, unitNormals1_y, unitNormals1_z;
    float areas1_i;
    float vector_P1P2_x, vector_P1P2_y, vector_P1P2_z, len_P1P2;
    float cos_alpha, cos_beta;
    float limbdarken, factor_single, gravitydarken;

    float g1_a0 = 1.0 / (g_ratio1 - 1.0 + tol);
    tmp1 = c1*c1 / (g_ratio2 - 1.0 + tol);
    tmp2 = a1*a1 / (g_ratio3 - 1.0 + tol);
    tmp3 = a1*a1*c1*c1;
    float g1_a2 = (tmp1 - tmp2) / (tmp3*(a1 + c1)) + g1_a0 * (a1 - c1) / tmp3;
    float g1_a1 = (tmp1 - g1_a0*c1*c1 - g1_a2*a1*tmp3) / tmp3;

    for(i=0; i<p1len; i++){
        if(thetas1[i] >= max_theta1)
            continue;
        i_ind2 = i*3+2;
        points1_z = points1[i_ind2];
        if(points1_z < 0.0)
            continue;
        i_ind0 = i*3; i_ind1 = i_ind0+1;
        points1_x = points1[i_ind0]; points1_y = points1[i_ind1]; 
        unitNormals1_x = unitNormals1[i_ind0]; unitNormals1_y = unitNormals1[i_ind1]; unitNormals1_z = unitNormals1[i_ind2];
        areas1_i = areas1[i];
        
        flag = 1;

        for(j=0; j<theta2len; j++){
            j_ind0 = j*2; j_ind1 = j_ind0 + 1;
            vector_P1P2_x = P2[j_ind0] - points1_x + sma;
            vector_P1P2_y = P2[j_ind1] - points1_y;
            vector_P1P2_z = - points1_z;

            tmp1 = vector_P1P2_x * unitNormals1_x;
            tmp2 = vector_P1P2_y * unitNormals1_y;
            tmp3 = vector_P1P2_z * unitNormals1_z;
            tmp1 += tmp2;
            cos_alpha = tmp1 + tmp3;
            if(cos_alpha <= 0.0)
                continue;

            tmp1 = vector_P1P2_x * np2[j_ind0];
            tmp2 = vector_P1P2_y * np2[j_ind1];
            cos_beta = tmp1 + tmp2;
            if(cos_beta >= 0.0)
                continue;

            if(flag){
                flag = 0;

                tmp1 = points1_y / b1;
                tmp1 *= tmp1;
                tmp2 = points1_x * points1_x;
                gravitydarken = 1.0 + (1.0 - tmp1) / (g1_a0 + g1_a1*tmp2 + g1_a2*tmp2*points1_x);
                gravitydarken = powf(gravitydarken, gdc1);
            }

            tmp1 = vector_P1P2_x * vector_P1P2_x;
            tmp2 = vector_P1P2_y * vector_P1P2_y;
            tmp3 = vector_P1P2_z * vector_P1P2_z;
            tmp1 += tmp2;
            tmp1 += tmp3;
            len_P1P2 = sqrtf(tmp1);
            cos_alpha /= len_P1P2;
            cos_beta /= len_P1P2;
            cos_beta *= -1.0;

            limbdarken = 1.0 - u1 * (1.0 - cos_alpha);
            factor_single = limbdarken * gravitydarken;
            factor_single *= cos_alpha;
            factor_single *= cos_beta;
            factor_single *= areas1_i;
            factor_single /= len_P1P2;
            factor_single /= len_P1P2;
            factor[j] += factor_single;
        }
    }
}

EXPORT void calc_np2_P2(const float* theta2, uint32_t theta2_len, float a, float b, float c, float* P2, float* np2){
    const float PI = 3.141592653589793;
    const float tol = 1.0e-6;

    uint32_t i;
    float theta2_single;
    float tmp1, tmp2, tmp3;

    float aa = a*a;
    float bb = b*b;
    float cc = c*c;
    float aa_bb = aa / bb;
    float cc_bb = cc / bb;
    float ab = a*b;
    float cb = c*b;
    uint32_t i_ind0, i_ind1;

    for(i=0; i<theta2_len; i++){
        theta2_single = theta2[i];
        i_ind0 = 2 * i;
        i_ind1 = i_ind0 + 1;
        if(theta2_single < PI/2.0 - tol){
            tmp1 = tanf(theta2_single);
            tmp2 = tmp1 * aa_bb;
            tmp3 = sqrtf(1.0 + tmp2 * tmp2);
            np2[i_ind0] = -1.0 / tmp3;
            np2[i_ind1] = tmp2 / tmp3;
            tmp2 = sqrtf(bb + aa * tmp1 * tmp1);
            P2[i_ind0] = - ab / tmp2;
            P2[i_ind1] = ab * tmp1 / tmp2;
        }
        else if(theta2_single > PI/2.0 + tol){
            tmp1 = tanf(theta2_single);
            tmp2 = tmp1 * cc_bb;
            tmp3 = sqrtf(1.0 + tmp2 * tmp2);
            np2[i_ind0] = 1.0 / tmp3;
            np2[i_ind1] = - tmp2 / tmp3;
            tmp2 = sqrtf(bb + cc * tmp1 * tmp1);
            P2[i_ind0] =  cb / tmp2;
            P2[i_ind1] = - cb * tmp1 / tmp2;
        }
        else{
            P2[i_ind0] = 0.0;
            P2[i_ind1] = b;
            np2[i_ind0] = 0.0;
            np2[i_ind1] = 1.0;
        }
    }
}


