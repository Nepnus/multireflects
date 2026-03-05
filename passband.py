import numpy as np
from scipy.interpolate import CubicSpline
import ctypes
from .clibs import lib_CubicSplineInterp as lib

class passband:
    u1 = None
    u2 = None
    A1 = None
    A2 = None
    
    bfac_interpfunc1 = None
    bfac_interpfunc2 = None
    a1_interpT = None
    b1_interpT = None
    c1_interpT = None
    d1_interpT = None
    a2_interpT = None
    b2_interpT = None
    c2_interpT = None
    d2_interpT = None
    h_logT = None
    logT_0 = None
    onestarflag = None

    def __init__(self, T, wave, I_lambda1, S_lambda, I_lambda2=None, T_interpnum=100):
        tol = 1.0e-4
        
        if I_lambda2 is None:
            self.onestarflag = True
        else:
            self.onestarflag = False
        
        S_lambda_brd = S_lambda[:, np.newaxis]
        S_lambda_brd = np.broadcast_to(S_lambda_brd, (len(S_lambda), len(T)))
        S_lambda_brd = np.transpose(S_lambda_brd, (1,0))
        IS_lambda1 = I_lambda1*S_lambda_brd
        flux1 = np.trapz(IS_lambda1, x=wave, axis=1)
        tmp_logT1 = 4.0*np.log(T)
        tmp_interpfunc = CubicSpline(tmp_logT1, flux1)
        tmp_logT2 = np.linspace(0.0, 1.0, T_interpnum, dtype=np.float32)
        tmp_logT2 /= tmp_logT2[-1]
        tmp_logT2 *= (tmp_logT1[-1] - tmp_logT1[0] - tol)
        tmp_logT2 += tmp_logT1[0]
        self.h_logT = tmp_logT2[1] - tmp_logT2[0]
        self.logT_0 = tmp_logT2[0]
        self.a1_interpT = np.zeros(T_interpnum, dtype=np.float32)
        y_interpT = tmp_interpfunc(tmp_logT2)
        np.copyto(self.a1_interpT, y_interpT)
        self.b1_interpT = np.zeros(T_interpnum-1, dtype=np.float32)
        self.c1_interpT = np.zeros(T_interpnum-1, dtype=np.float32)
        self.d1_interpT = np.zeros(T_interpnum-1, dtype=np.float32)
        lib.getInterpFactors(
            self.a1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            self.h_logT, len(tmp_logT2),
            self.b1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            self.c1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            self.d1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        )
        
        if self.onestarflag == False:
            IS_lambda2 = I_lambda2*S_lambda_brd
            flux2 = np.trapz(IS_lambda2, x=wave, axis=1)
            tmp_interpfunc = CubicSpline(tmp_logT1, flux2)
            self.a2_interpT = np.zeros(T_interpnum, dtype=np.float32)
            y_interpT = tmp_interpfunc(tmp_logT2)
            np.copyto(self.a2_interpT, y_interpT)
            self.b2_interpT = np.zeros(T_interpnum-1, dtype=np.float32)
            self.c2_interpT = np.zeros(T_interpnum-1, dtype=np.float32)
            self.d2_interpT = np.zeros(T_interpnum-1, dtype=np.float32)
            lib.getInterpFactors(
                self.a2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.h_logT, len(tmp_logT2),
                self.b2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.c2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.d2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )
        
        wave_brd = wave[:, np.newaxis]
        wave_brd = np.broadcast_to(wave_brd, (len(wave), len(T)))
        wave_brd = np.transpose(wave_brd, (1,0))
        diff_I_lambda1 = np.diff(I_lambda1, axis=1)
        diff_wave = np.diff(wave_brd, axis=1)
        bfac_lambda1 = 5.0 + wave_brd[:,1:] / I_lambda1[:,1:] * diff_I_lambda1 / diff_wave
        bfac_list1 = np.trapz(IS_lambda1[:,1:]*bfac_lambda1, x=wave[1:], axis=1) / np.trapz(IS_lambda1[:,1:], x=wave[1:], axis=1)
        ind1 = ~np.isnan(bfac_list1)
        self.bfac_interpfunc1 = CubicSpline(T[ind1], bfac_list1[ind1])
        
        if self.onestarflag == False:
            diff_I_lambda2 = np.diff(I_lambda2, axis=1)
            bfac_lambda2 = 5.0 + wave_brd[:,1:] / I_lambda2[:,1:] * diff_I_lambda2 / diff_wave
            bfac_list2 = np.trapz(IS_lambda2[:,1:]*bfac_lambda2, x=wave[1:], axis=1) / np.trapz(IS_lambda2[:,1:], x=wave[1:], axis=1)
            ind2 = ~np.isnan(bfac_list2)
            self.bfac_interpfunc2 = CubicSpline(T[ind2], bfac_list2[ind2])
    
    def get_bfac(self, T1, T2=None):
        if self.onestarflag:
            return self.bfac_interpfunc1(T1)
        else:
            return self.bfac_interpfunc1(T1), self.bfac_interpfunc2(T2)


