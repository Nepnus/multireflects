#include "CubicSplineInterp.h"
#include <math.h>

void getInterpFactors(const fn* y, fn h, uint32_t len, fn* b, fn* c, fn* d){
    const fn sqrt_3 = (fn)sqrt(3.0);

    fn m, n;
    fn f, g, s2, s1;
    uint32_t i;
    fn power_tmp = 1.0;

    f = 6.0 / h / h * (y[2] - 2.0*y[1] + y[0]);
    c[1] = f * 0.25;

    if(len>=4)
        for(i=2; i<len-1; i++){
            f = 6.0 / h / h * (y[i+1] - 2.0*y[i] + y[i-1]);
            if(i <= 15){
                n = (97.0 + 56.0*sqrt_3) * power_tmp;
                m = (n * (2.0 - sqrt_3) - 2.0 - sqrt_3) / (n - 1.0);
                power_tmp *= (7.0 + 4.0*sqrt_3);
            }
            c[i] = (f - c[i-1]) / (4.0 - m);
        }

    s1 = c[len-2];
    b[len-2] = (y[len-1] - y[len-2]) / h - h * s1 / 3.0;
    d[len-2] = -s1 / h / 6.0;
    c[len-2] /= 2.0;

    if(len>=4)
        for(i=len-3; i>=1; i--){
            if(i <= 15){
                power_tmp /= (7.0 + 4.0*sqrt_3);
                n = (97.0 + 56.0*sqrt_3) * power_tmp;
                m = (n * (2.0 - sqrt_3) - 2.0 - sqrt_3) / (n - 1.0);
            }
            s2 = c[i+1] * 2.0;
            s1 = c[i] - m * s2;
            c[i] = s1 / 2.0;
            b[i] = (y[i+1] - y[i]) / h - h * (0.5 * s1 + (s2 - s1) / 6.0);
            d[i] = (s2 - s1) / h / 6.0;
        }

    s2 = c[1] * 2.0;
    b[0] = (y[1] - y[0]) / h - h * s2 / 6.0;
    c[0] = 0.0;
    d[0] = s2 / h / 6.0;
}

void getInterpFactors_bycol(const fn* y, fn h, uint32_t len, uint32_t rowlen, fn* a, fn* b, fn* c, fn* d){
    const fn sqrt_3 = (fn)sqrt(3.0);

    fn m, n;
    fn f, g, s2, s1;
    uint32_t i;
    fn power_tmp = 1.0;

    f = 6.0 / h / h * (y[2*rowlen] - 2.0*y[rowlen] + y[0]);
    c[1] = f * 0.25;

    if(len>=4)
        for(i=2; i<len-1; i++){
            f = 6.0 / h / h * (y[(i+1)*rowlen] - 2.0*y[i*rowlen] + y[(i-1)*rowlen]);
            if(i <= 15){
                n = (97.0 + 56.0*sqrt_3) * power_tmp;
                m = (n * (2.0 - sqrt_3) - 2.0 - sqrt_3) / (n - 1.0);
                power_tmp *= (7.0 + 4.0*sqrt_3);
            }
            c[i] = (f - c[i-1]) / (4.0 - m);
        }

    s1 = c[len-2];
    b[len-2] = (y[(len-1)*rowlen] - y[(len-2)*rowlen]) / h - h * s1 / 3.0;
    d[len-2] = -s1 / h / 6.0;
    c[len-2] /= 2.0;
    a[len-2] = y[(len-2)*rowlen];

    if(len>=4)
        for(i=len-3; i>=1; i--){
            if(i <= 15){
                power_tmp /= (7.0 + 4.0*sqrt_3);
                n = (97.0 + 56.0*sqrt_3) * power_tmp;
                m = (n * (2.0 - sqrt_3) - 2.0 - sqrt_3) / (n - 1.0);
            }
            s2 = c[i+1] * 2.0;
            s1 = c[i] - m * s2;
            c[i] = s1 / 2.0;
            b[i] = (y[(i+1)*rowlen] - y[i*rowlen]) / h - h * (0.5 * s1 + (s2 - s1) / 6.0);
            d[i] = (s2 - s1) / h / 6.0;
            a[i] = y[i*rowlen];
        }

    s2 = c[1] * 2.0;
    b[0] = (y[rowlen] - y[0]) / h - h * s2 / 6.0;
    c[0] = 0.0;
    d[0] = s2 / h / 6.0;
    a[0] = y[0];
}

fn interp_single(fn x_interp, fn h, fn x_0, const fn* a, const fn* b, const fn* c, const fn* d){
    uint32_t a_i;
    fn y_interp, tmp1;

    a_i = (uint32_t) ((x_interp - x_0) / h);
    tmp1 = x_interp - a_i*h - x_0;
    y_interp = a[a_i] + b[a_i]*tmp1 + c[a_i]*tmp1*tmp1 + d[a_i]*tmp1*tmp1*tmp1;
    return y_interp;
}

void interp(const fn* x_interp, uint32_t len, fn h, fn x_0, const fn* a, const fn* b, const fn* c, const fn* d, fn* y_interp){
    uint32_t x_i, a_i;
    fn tmp1;

    for(x_i=0; x_i<len; x_i++){
        a_i = (uint32_t) ((x_interp[x_i] - x_0) / h);
        tmp1 = x_interp[x_i] - a_i*h - x_0;
        y_interp[x_i] = a[a_i] + b[a_i]*tmp1 + c[a_i]*tmp1*tmp1 + d[a_i]*tmp1*tmp1*tmp1;
    }

}

