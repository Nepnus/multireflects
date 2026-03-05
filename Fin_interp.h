#include <stdint.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#ifndef __Fin_interp__
#define __Fin_interp__

EXPORT void get_ddirection_factors(const float* factors, uint32_t theta2_num, uint32_t d_num, float d_h, float* A, float* B, float* C, float* D);
EXPORT void ddirection_interp(float d_interp, float d_h, float d_0, uint32_t theta2_num, uint32_t d_num,
    const float* A, const float* B, const float* C, const float* D, float* y);
EXPORT void get_thetadirection_factors(const float* factors, float theta_h, uint32_t theta2_num, float* b, float* c, float* d);
EXPORT float thetadirection_interp(float theta_interp, float theta_h, float theta_max,
    const float* a, const float* b, const float* c, const float* d);

#endif

