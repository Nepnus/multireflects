import numpy as np
from scipy.spatial import ConvexHull
import ctypes
from .clibs import lib_UnitsCore as lib

class surfaceUnits:
    thetaInterval = None
    phiInterval = None
    a = None
    b = None
    b2 = None
    c = None
    faces = None
    points_surface = None
    
    centroids = None
    areas = None
    unitNormals = None
    thetas = None
    
    def __init__(self, thetaInterval=0.05, phiInterval=0.05):
        theta = np.arange(0.0, np.pi+thetaInterval, thetaInterval, dtype=np.float32)
        ind = (theta>=0.0) & (theta<=np.pi)
        theta = theta[ind]
        for i in range(len(theta)):
            theta_single = theta[i]
            phiInterval_theta = phiInterval / np.cos(np.pi/2.0 - theta_single)
            phi_theta = np.arange(0.0, 2.0*np.pi+phiInterval_theta, phiInterval_theta, dtype=np.float32)
            ind = (phi_theta>=0) & (phi_theta<2.0*np.pi)
            phi_theta = phi_theta[ind]
            
            x = np.sin(theta_single) * np.sin(phi_theta)
            y = np.cos(theta_single) * np.ones_like(x, dtype=np.float32)
            z = np.sin(theta_single) * np.cos(phi_theta)
            points_theta = np.array([x,y,z], dtype=np.float32)
            points_theta = points_theta.T
            
            if self.points_surface is None:
                self.points_surface = points_theta
            else:
                self.points_surface = np.concatenate([self.points_surface, points_theta], axis=0)
        
        self.points_surface = self.points_surface.copy()
        hull = ConvexHull(self.points_surface)
        self.faces = hull.simplices

        self.centroids = np.zeros_like(self.faces, dtype=np.float32)
        self.areas = np.zeros(self.faces.shape[0], dtype=np.float32)
        self.unitNormals = np.zeros_like(self.faces, dtype=np.float32)
        self.thetas = np.zeros(self.faces.shape[0], dtype=np.float32)
    
    def get_areaUnits(self, a, b, b2, c, thetasflag=True):
        self.a = a; self.b = b; self.b2 = b2; self.c = c
        if thetasflag:
            lib.core_thetas(
                self.points_surface.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.faces.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                self.faces.shape[0], a, b, b2, c,
                self.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.areas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.unitNormals.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.thetas.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )
        else:
            lib.core_nothetas(
                self.points_surface.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.faces.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                self.faces.shape[0], a, b, b2, c,
                self.centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.areas.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                self.unitNormals.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            )

