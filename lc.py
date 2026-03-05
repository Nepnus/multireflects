import numpy as np
import ctypes
from . import calc_distance
from . import calc_shape
from . import get_areaUnits
from . import calc_Fin
from .clibs import lib_lcmode0 as lib0
from .clibs import lib_lcmode1 as lib1

class lc_workflow:
    tol = 1.0e-4
    shape1 = None; shape2 = None
    T1 = None; T2 = None
    r1 = None; r2 = None
    A1 = None; A2 = None
    gdc1 = None; gdc2 = None
    u1 = None; u2 = None
    Fin1 = None; Fin2 = None
    
    phi = None; nv = None
    incl = None; q = None
    P = None; sma = None
    ecc = None; w2 = None
    RV1 = None; RV2 = None
    RV_gamma = None
    d = None

    lc = None
    lc1 = None; lc2 = None

    cos_gamma1 = None; cos_gamma2 = None
    Teff_distribution1 = None; Teff_distribution2 = None
    incident_light_factor1 = None; incident_light_factor2 = None
    tmp1 = None; tmp2 = None
    i_ps1 = None; i_sum1 = None
    i_ps2 = None; i_sum2 = None
    
    def __init__(self, phi, passbandnum=1, thetaInterval=0.05, phiInterval=0.05, theta2_num=10, d_num=5):
        self.phi = phi.copy()
        self.shape1 = get_areaUnits.surfaceUnits(thetaInterval, phiInterval)
        self.shape2 = get_areaUnits.surfaceUnits(thetaInterval, phiInterval)
        self.Fin1 = calc_Fin.Fin_workflow(self.shape2, self.shape1, theta2_num, d_num)
        self.Fin2 = calc_Fin.Fin_workflow(self.shape1, self.shape2, theta2_num, d_num)
        self.lc = np.zeros((passbandnum, len(phi)), dtype=np.float32)
        self.lc1 = np.zeros((passbandnum, len(phi)), dtype=np.float32)
        self.lc2 = np.zeros((passbandnum, len(phi)), dtype=np.float32)
        
        self.cos_gamma1 = np.zeros(len(phi)*self.shape1.centroids.shape[0], dtype=np.float32)
        self.cos_gamma2 = np.zeros_like(self.cos_gamma1, dtype=np.float32)
        self.Teff_distribution1 = np.zeros_like(self.cos_gamma1, dtype=np.float32)
        self.Teff_distribution2 = np.zeros_like(self.cos_gamma1, dtype=np.float32)
        self.incident_light_factor1 = np.zeros_like(self.cos_gamma1, dtype=np.float32)
        self.incident_light_factor2 = np.zeros_like(self.cos_gamma2, dtype=np.float32)
        self.tmp1 = np.zeros(self.shape1.centroids.shape[0], dtype=np.float32)
        self.tmp2 = np.zeros(self.shape1.centroids.shape[0], dtype=np.float32)
        self.i_ps1 = np.zeros_like(self.cos_gamma1, dtype=np.uint32)
        self.i_ps2 = np.zeros_like(self.cos_gamma1, dtype=np.uint32)
        self.i_sum1 = np.zeros(len(phi), dtype=np.uint32)
        self.i_sum2 = np.zeros(len(phi), dtype=np.uint32)

    def lc_mode1(self, passbandlist):
        self.d, self.nv = calc_distance.calc_d(self.phi, self.ecc, self.w2)
        self.d *= self.sma
        self.RV1, self.RV2 = calc_distance.calc_rv(self.nv, self.ecc, self.w2, self.incl, self.q, self.sma, self.P)
        self.RV1 += self.RV_gamma
        self.RV2 += self.RV_gamma
        
        a1, b12, c1, g_abratio1, g_b2bratio1, g_cbratio1 = calc_shape.get_shape(1.0/self.q, self.sma, self.r1)
        a2, b22, c2, g_abratio2, g_b2bratio2, g_cbratio2 = calc_shape.get_shape(self.q, self.sma, self.r2)
        self.shape1.get_areaUnits(c1, self.r1, b12, a1, False)
        self.shape2.get_areaUnits(c2, self.r2, b22, a2)
        
        self.Fin2.get_F(self.sma, self.ecc, self.u1, self.gdc1, g_abratio1, g_b2bratio1, g_cbratio1)
        
        lib0.gen_cosgamma_multid(
            self.incl, self.nv.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            len(self.nv), a1, self.r1, b12, c1, a2, self.r2, b22, c2, 
            self.shape1.unitNormals.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.unitNormals.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
            self.shape1.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            self.shape1.centroids.shape[0], self.shape2.centroids.shape[0], 
            self.cos_gamma1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.i_ps1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)), self.i_sum1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)), 
            self.cos_gamma2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.i_ps2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)), self.i_sum2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        )
        
        lib1.gen_Teff_distribution(
            self.shape1.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape1.centroids.shape[0],
            self.T1, self.gdc1, 
            self.shape1.c, self.shape1.b, self.shape1.a,
            g_abratio1, g_b2bratio1, g_cbratio1,
            self.Teff_distribution1.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        )
        
        if self.ecc < self.tol:
            lib0.gen_Teff_distribution_singled(
                self.shape2.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.thetas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.centroids.shape[0],
                self.T2, self.T1, self.A2, self.gdc2,
                self.shape2.c, self.shape2.b, self.shape2.a,
                g_abratio2, g_b2bratio2, g_cbratio2,
                self.Fin2.theta_h, self.Fin2.max_theta2, 
                self.Fin2.factor.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin2.c.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.di.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Teff_distribution2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor2.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )
            
            for i in range(len(passbandlist)):
                passband_obj = passbandlist[i]
                lc1_row = self.lc1[i, :]; lc2_row = self.lc2[i, :]
                lc1_row_ptr = ctypes.cast(lc1_row.ctypes.data, ctypes.POINTER(ctypes.c_float))
                lc2_row_ptr = ctypes.cast(lc2_row.ctypes.data, ctypes.POINTER(ctypes.c_float))
                bfac1, bfac2 = passband_obj.get_bfac(self.T1, self.T2)
                
                lib1.gen_lc(
                    self.Teff_distribution1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.u1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.shape1.areas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.cos_gamma1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.RV1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), bfac1, 
                    passband_obj.a1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.h_logT, passband_obj.logT_0, 
                    self.i_ps1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    len(self.phi), lc1_row_ptr, self.tmp1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(self.shape1.areas)
                )
                lib0.gen_lc_singled(
                    self.Teff_distribution2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    passband_obj.u2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.areas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    self.cos_gamma2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.RV2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), bfac2, passband_obj.A2, self.T1,
                    passband_obj.a2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.h_logT, passband_obj.logT_0, 
                    passband_obj.a1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    self.i_ps2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    len(self.phi), lc2_row_ptr, self.tmp1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(self.shape2.areas)
                )
                np.add(lc1_row, lc2_row, out=self.lc[i, :])
        else:
            lib0.gen_Teff_distribution_multid(
                self.shape2.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.thetas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.centroids.shape[0],
                self.T2, self.T1, self.A2, self.gdc2,
                self.shape2.c, self.shape2.b, self.shape2.a,
                g_abratio2, g_b2bratio2, g_cbratio2,
                self.d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(self.d), self.i_ps2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                self.Fin2.d_h, self.Fin2.d[0], self.Fin2.theta_h, self.Fin2.max_theta2, self.Fin2.d_num, self.Fin2.theta2_num,
                self.Fin2.A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.B.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin2.C.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.D.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin2.a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin2.c.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.di.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Teff_distribution2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.tmp1.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )
            
            for i in range(len(passbandlist)):
                passband_obj = passbandlist[i]
                lc1_row = self.lc1[i, :]; lc2_row = self.lc2[i, :]
                lc1_row_ptr = ctypes.cast(lc1_row.ctypes.data, ctypes.POINTER(ctypes.c_float))
                lc2_row_ptr = ctypes.cast(lc2_row.ctypes.data, ctypes.POINTER(ctypes.c_float))
                bfac1, bfac2 = passband_obj.get_bfac(self.T1, self.T2)
                
                lib1.gen_lc(
                    self.Teff_distribution1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.u1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.shape1.areas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.cos_gamma1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.RV1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), bfac1, 
                    passband_obj.a1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.h_logT, passband_obj.logT_0, 
                    self.i_ps1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    len(self.phi), lc1_row_ptr, self.tmp1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(self.shape1.areas)
                )
                lib0.gen_lc_multid(
                    self.Teff_distribution2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    passband_obj.u2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.areas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    self.cos_gamma2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.RV2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), bfac2, passband_obj.A2, self.T1,
                    passband_obj.a2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.h_logT, passband_obj.logT_0, 
                    passband_obj.a1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    self.i_ps2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    len(self.phi), lc2_row_ptr
                )
                np.add(lc1_row, lc2_row, out=self.lc[i, :])
        
    def lc_mode0(self, passbandlist):
        self.d, self.nv = calc_distance.calc_d(self.phi, self.ecc, self.w2)
        self.d *= self.sma
        self.RV1, self.RV2 = calc_distance.calc_rv(self.nv, self.ecc, self.w2, self.incl, self.q, self.sma, self.P)
        self.RV1 += self.RV_gamma
        self.RV2 += self.RV_gamma
        
        a1, b12, c1, g_abratio1, g_b2bratio1, g_cbratio1 = calc_shape.get_shape(1.0/self.q, self.sma, self.r1)
        a2, b22, c2, g_abratio2, g_b2bratio2, g_cbratio2 = calc_shape.get_shape(self.q, self.sma, self.r2)
        self.shape1.get_areaUnits(c1, self.r1, b12, a1)
        self.shape2.get_areaUnits(c2, self.r2, b22, a2)
        
        self.Fin1.get_F(self.sma, self.ecc, self.u2, self.gdc2, g_abratio2, g_b2bratio2, g_cbratio2)
        self.Fin2.get_F(self.sma, self.ecc, self.u1, self.gdc1, g_abratio1, g_b2bratio1, g_cbratio1)
        
        lib0.gen_cosgamma_multid(
            self.incl, self.nv.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            len(self.nv), a1, self.r1, b12, c1, a2, self.r2, b22, c2, 
            self.shape1.unitNormals.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.unitNormals.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
            self.shape1.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            self.shape1.centroids.shape[0], self.shape2.centroids.shape[0], 
            self.cos_gamma1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.i_ps1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)), self.i_sum1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)), 
            self.cos_gamma2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.i_ps2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)), self.i_sum2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        )
        
        if self.ecc < self.tol:
            lib0.gen_Teff_distribution_singled(
                self.shape1.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape1.thetas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape1.centroids.shape[0],
                self.T1, self.T2, self.A1, self.gdc1,
                self.shape1.c, self.shape1.b, self.shape1.a, 
                g_abratio1, g_b2bratio1, g_cbratio1,
                self.Fin1.theta_h, self.Fin1.max_theta2, 
                self.Fin1.factor.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin1.b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin1.c.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin1.di.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Teff_distribution1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor1.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )
            lib0.gen_Teff_distribution_singled(
                self.shape2.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.thetas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.centroids.shape[0],
                self.T2, self.T1, self.A2, self.gdc2,
                self.shape2.c, self.shape2.b, self.shape2.a, 
                g_abratio2, g_b2bratio2, g_cbratio2,
                self.Fin2.theta_h, self.Fin2.max_theta2, 
                self.Fin2.factor.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin2.c.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.di.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Teff_distribution2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor2.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )
            
            for i in range(len(passbandlist)):
                passband_obj = passbandlist[i]
                lc1_row = self.lc1[i, :]; lc2_row = self.lc2[i, :]
                lc1_row_ptr = ctypes.cast(lc1_row.ctypes.data, ctypes.POINTER(ctypes.c_float))
                lc2_row_ptr = ctypes.cast(lc2_row.ctypes.data, ctypes.POINTER(ctypes.c_float))
                bfac1, bfac2 = passband_obj.get_bfac(self.T1, self.T2)
                
                lib0.gen_lc_singled(
                    self.Teff_distribution1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    passband_obj.u1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape1.areas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.cos_gamma1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.RV1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), bfac1, passband_obj.A1, self.T2,
                    passband_obj.a1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.h_logT, passband_obj.logT_0, 
                    passband_obj.a2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.i_ps1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    len(self.phi), lc1_row_ptr, self.tmp1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(self.shape1.areas)
                )
                lib0.gen_lc_singled(
                    self.Teff_distribution2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    passband_obj.u2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.areas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.cos_gamma2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.RV2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), bfac2, passband_obj.A2, self.T1,
                    passband_obj.a2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.h_logT, passband_obj.logT_0, 
                    passband_obj.a1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.i_ps2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    len(self.phi), lc2_row_ptr, self.tmp1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(self.shape2.areas)
                )
                np.add(lc1_row, lc2_row, out=self.lc[i, :])
        else:
            lib0.gen_Teff_distribution_multid(
                self.shape1.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape1.thetas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape1.centroids.shape[0],
                self.T1, self.T2, self.A1, self.gdc1,
                self.shape1.c, self.shape1.b, self.shape1.a, 
                g_abratio1, g_b2bratio1, g_cbratio1,
                self.d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(self.d), self.i_ps1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                self.Fin1.d_h, self.Fin1.d[0], self.Fin1.theta_h, self.Fin1.max_theta2, self.Fin1.d_num, self.Fin1.theta2_num,
                self.Fin1.A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin1.B.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin1.C.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin1.D.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin1.a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin1.b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin1.c.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin1.di.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Teff_distribution1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.tmp1.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )
            lib0.gen_Teff_distribution_multid(
                self.shape2.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.thetas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.centroids.shape[0],
                self.T2, self.T1, self.A2, self.gdc2,
                self.shape2.c, self.shape2.b, self.shape2.a, 
                g_abratio2, g_b2bratio2, g_cbratio2,
                self.d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(self.d), self.i_ps2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                self.Fin2.d_h, self.Fin2.d[0], self.Fin2.theta_h, self.Fin2.max_theta2, self.Fin2.d_num, self.Fin2.theta2_num,
                self.Fin2.A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.B.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin2.C.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.D.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin2.a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Fin2.c.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.Fin2.di.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                self.Teff_distribution2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.tmp1.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )

            for i in range(len(passbandlist)):
                passband_obj = passbandlist[i]
                lc1_row = self.lc1[i, :]; lc2_row = self.lc2[i, :]
                lc1_row_ptr = ctypes.cast(lc1_row.ctypes.data, ctypes.POINTER(ctypes.c_float))
                lc2_row_ptr = ctypes.cast(lc2_row.ctypes.data, ctypes.POINTER(ctypes.c_float))
                bfac1, bfac2 = passband_obj.get_bfac(self.T1, self.T2)

                lib0.gen_lc_multid(
                    self.Teff_distribution1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    passband_obj.u1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape1.areas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.cos_gamma1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.RV1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), bfac1, passband_obj.A1, self.T2,
                    passband_obj.a1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.h_logT, passband_obj.logT_0, 
                    passband_obj.a2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.i_ps1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum1.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    len(self.phi), lc1_row_ptr
                )
                lib0.gen_lc_multid(
                    self.Teff_distribution2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.incident_light_factor2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    passband_obj.u2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), self.shape2.areas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    self.cos_gamma2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.RV2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), bfac2,  passband_obj.A2, self.T1,
                    passband_obj.a2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d2_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.h_logT, passband_obj.logT_0, 
                    passband_obj.a1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.b1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), 
                    passband_obj.c1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), passband_obj.d1_interpT.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    self.i_ps2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),  self.i_sum2.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                    len(self.phi), lc2_row_ptr
                )
                np.add(lc1_row, lc2_row, out=self.lc[i, :])

