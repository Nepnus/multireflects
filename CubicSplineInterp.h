#include <stdint.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

#ifndef __CubicSplineInterp__
#define __CubicSplineInterp__

typedef float fn;

EXPORT void getInterpFactors(const fn* y, fn h, uint32_t len, fn* b, fn* c, fn* d);
EXPORT void interp(const fn* x_interp, uint32_t len, fn h, fn x_0, const fn* a, const fn* b, const fn* c, const fn* d, fn* y_interp);
EXPORT void getInterpFactors_bycol(const fn* y, fn h, uint32_t len, uint32_t rowlen, fn* a, fn* b, fn* c, fn* d);
EXPORT fn interp_single(fn x_interp, fn h, fn x_0, const fn* a, const fn* b, const fn* c, const fn* d);

#endif
