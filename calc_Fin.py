import numpy as np
import ctypes
from .clibs import lib_FinCore as lib
from .clibs import lib_FinInterp as lib2

class Fin_workflow:
    tol = 1.0e-4
    max_theta2 = 0.0
    ind = None
    ind_num = 0
    shape1 = None
    shape2 = None
    
    theta2_num = 10; d_num = 5
    theta2 = None; d = None
    d_singleflag = False
    factor = None
    interpfunc = None
    
    d_h = 0.0; theta_h = 0.0
    a = None
    b = None
    c = None
    di = None
    A = None
    B = None
    C = None
    D = None
    
    P2 = None
    n_p2 = None
    
    def __init__(self, shape1, shape2, theta2_num=10, d_num=5):
        self.theta2_num = theta2_num; self.d_num = d_num
        self.theta2 = np.linspace(0.0, 1.0, theta2_num, dtype=np.float32)
        self.theta2 /= self.theta2[-1]
        self.d = np.linspace(0.0, 1.0, d_num, dtype=np.float32)
        self.d /= self.d[-1]
        self.shape1 = shape1
        self.shape2 = shape2
        self.factor = np.zeros((d_num, theta2_num), dtype=np.float32)
        self.a = np.zeros(theta2_num, dtype=np.float32)
        self.b = np.zeros(theta2_num-1, dtype=np.float32)
        self.c = np.zeros(theta2_num-1, dtype=np.float32)
        self.di = np.zeros(theta2_num-1, dtype=np.float32)
        self.A = np.zeros((theta2_num, d_num-1), dtype=np.float32)
        self.B = np.zeros((theta2_num, d_num-1), dtype=np.float32)
        self.C = np.zeros((theta2_num, d_num-1), dtype=np.float32)
        self.D = np.zeros((theta2_num, d_num-1), dtype=np.float32)
        
        self.P2 = np.zeros((theta2_num, 2), dtype=np.float32)
        self.n_p2 = np.zeros((theta2_num, 2), dtype=np.float32)
    
    def get_F(self, sma, ecc, u1, gdc1, g_ratio1, g_ratio2, g_ratio3):
        if ecc < self.tol:
            self.get_F_single_d(sma, u1, gdc1, g_ratio1, g_ratio2, g_ratio3)
            self.d_singleflag = True
            lib2.get_thetadirection_factors(
                self.factor.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.theta_h, self.theta2_num,
                self.b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.c.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.di.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )
        else:
            self.get_F_multi_d(sma, ecc, u1, gdc1, g_ratio1, g_ratio2, g_ratio3)
            self.d_singleflag = False
            lib2.get_ddirection_factors(
                self.factor.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.theta2_num, self.d_num, self.d_h,
                self.A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.B.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.C.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.D.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )
    
    def get_F_single_d(self, sma, u1, gdc1, g_ratio1, g_ratio2, g_ratio3):
        shape1 = self.shape1; shape2 = self.shape2
        max_theta2 = np.arctan2(np.sqrt(sma**2.0 - (shape2.b-shape1.b)**2.0), (shape2.b-shape1.b))
        self.max_theta2 = max_theta2
        max_theta1 = np.pi - max_theta2

        self.theta2 /= self.theta2[-1]
        self.theta2 *= (max_theta2 + self.tol)
        self.theta_h = self.theta2[1] - self.theta2[0]
        self.factor[0, :] = 0.0
        
        points1 = shape1.centroids
        areas1 = shape1.areas
        unitNormals1 = shape1.unitNormals
        
        theta2 = self.theta2
        P2 = self.P2
        n_p2 = self.n_p2
        lib.calc_np2_P2(
            theta2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            len(theta2),
            shape2.c, shape2.b, shape2.a,
            P2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            n_p2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        )
        
        lib.matrixcore(
            points1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            areas1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            unitNormals1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            shape1.thetas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            P2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            n_p2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            points1.shape[0], len(theta2),
            max_theta1, sma, u1, gdc1, shape1.c, shape1.b, shape1.a,
            g_ratio1, g_ratio2, g_ratio3,
            self.factor.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        )
        self.factor[0, :] *= 2.0
        
    def get_F_multi_d(self, sma, ecc, u1, gdc1, g_ratio1, g_ratio2, g_ratio3):
        shape1 = self.shape1; shape2 = self.shape2
        max_theta2 = np.arctan2(np.sqrt(sma**2.0 - (shape2.b-shape1.b)**2.0), (shape2.b-shape1.b))
        self.max_theta2 = max_theta2
        max_theta1 = np.pi - max_theta2
        
        self.theta2 /= self.theta2[-1]
        self.theta2 *= (max_theta2 + self.tol)
        self.d -= self.d[0]
        self.d /= self.d[-1]
        self.d *= (2.0*sma*(ecc+self.tol))
        self.d += (sma*(1.0-ecc-self.tol))
        self.d_h = self.d[1] - self.d[0]
        self.theta_h = self.theta2[1] - self.theta2[0]
        self.factor[:, :] = 0.0
        
        points1 = shape1.centroids
        areas1 = shape1.areas
        unitNormals1 = shape1.unitNormals
        
        theta2 = self.theta2
        P2 = self.P2
        n_p2 = self.n_p2
        lib.calc_np2_P2(
            theta2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            len(theta2),
            shape2.c, shape2.b, shape2.a,
            P2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            n_p2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        )
        
        lib.tensorcore(
            points1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            areas1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            unitNormals1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            shape1.thetas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            P2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            n_p2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            self.d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            points1.shape[0], len(theta2), len(self.d),
            max_theta1, u1, gdc1, shape1.c, shape1.b, shape1.a,
            g_ratio1, g_ratio2, g_ratio3,
            self.factor.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        )
        self.factor[:, :] *= 2.0


