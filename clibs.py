import ctypes
import os
import sys

libs_suffix = None
if sys.platform.startswith("win"):
    libs_suffix = ".dll"
elif sys.platform.startswith("linux"):
    libs_suffix = ".so"
elif sys.platform.startswith("darwin"):
    libs_suffix = ".dylib"
else:
    raise Exception("Unsupportted system.")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

lib_CubicSplineInterp_path = "libCubicSplineInterp"
lib_CubicSplineInterp_path += libs_suffix
lib_CubicSplineInterp_path = os.path.join(BASE_DIR, lib_CubicSplineInterp_path)
lib_CubicSplineInterp = ctypes.CDLL(lib_CubicSplineInterp_path, mode=ctypes.RTLD_GLOBAL)
lib_CubicSplineInterp.getInterpFactors.argtypes = [
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_float, ctypes.c_uint32,
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float)
]
lib_CubicSplineInterp.getInterpFactors.restype = None
lib_CubicSplineInterp.interp_single.argtypes = [
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)
]
lib_CubicSplineInterp.interp_single.restype = ctypes.c_float

lib_FinInterp_path = "libFin_interp"
lib_FinInterp_path += libs_suffix
lib_FinInterp_path = os.path.join(BASE_DIR, lib_FinInterp_path)
lib_FinInterp = ctypes.CDLL(lib_FinInterp_path, mode=ctypes.RTLD_GLOBAL)
lib_FinInterp.get_ddirection_factors.argtypes = [
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_uint32, ctypes.c_uint32, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float)
]
lib_FinInterp.get_ddirection_factors.restype = None
lib_FinInterp.get_thetadirection_factors.argtypes = [
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_float, ctypes.c_uint32,
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float)
]
lib_FinInterp.get_thetadirection_factors.restype = None
    
lib_UnitsCore_path = "Units_core"
lib_UnitsCore_path += libs_suffix
lib_UnitsCore_path = os.path.join(BASE_DIR, lib_UnitsCore_path)
lib_UnitsCore = ctypes.CDLL(lib_UnitsCore_path)
lib_UnitsCore.core_nothetas.argtypes = [
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_uint32),
    ctypes.c_uint32,
    ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float)
]
lib_UnitsCore.core_nothetas.restype = None
lib_UnitsCore.core_thetas.argtypes = [
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_uint32),
    ctypes.c_uint32,
    ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float)
]
lib_UnitsCore.core_thetas.restype = None

lib_FinCore_path = "Fin_core"
lib_FinCore_path += libs_suffix
lib_FinCore_path = os.path.join(BASE_DIR, lib_FinCore_path)
lib_FinCore = ctypes.CDLL(lib_FinCore_path)
lib_FinCore.matrixcore.argtypes = [
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_uint32, ctypes.c_uint32,
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float)
]
lib_FinCore.matrixcore.restype = None
lib_FinCore.tensorcore.argtypes = [
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_uint32, ctypes.c_uint32, ctypes.c_uint32,
    ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float)
]
lib_FinCore.tensorcore.restype = None
lib_FinCore.calc_np2_P2.argtypes = [
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_uint32,
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float)
]
lib_FinCore.calc_np2_P2.restype = None

lib_lcmode0_path = "lc_mode0"
lib_lcmode0_path += libs_suffix
lib_lcmode0_path = os.path.join(BASE_DIR, lib_lcmode0_path)
lib_lcmode0 = ctypes.CDLL(lib_lcmode0_path)
lib_lcmode0.gen_cosgamma_multid.argtypes = [
    ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.c_uint32,
    ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.c_uint32, ctypes.c_uint32,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_uint32),
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_uint32)
]
lib_lcmode0.gen_cosgamma_multid.restype = None
lib_lcmode0.gen_Teff_distribution_multid.argtypes = [
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.c_uint32,
    ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float), ctypes.c_uint32, ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_uint32),
    ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_uint32, ctypes.c_uint32,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)
]
lib_lcmode0.gen_Teff_distribution_multid.restype = None
lib_lcmode0.gen_lc_multid.argtypes = [
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float), ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_uint32), ctypes.c_uint32, ctypes.POINTER(ctypes.c_float)
]
lib_lcmode0.gen_lc_multid.restype = None
lib_lcmode0.gen_Teff_distribution_singled.argtypes = [
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.c_uint32,
    ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)
]
lib_lcmode0.gen_Teff_distribution_singled.restype = None
lib_lcmode0.gen_lc_singled.argtypes = [
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float), ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_uint32), ctypes.c_uint32, ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float), ctypes.c_uint32
]
lib_lcmode0.gen_lc_singled.restype = None

lib_lcmode1_path = "lc_mode1"
lib_lcmode1_path += libs_suffix
lib_lcmode1_path = os.path.join(BASE_DIR, lib_lcmode1_path)
lib_lcmode1 = ctypes.CDLL(lib_lcmode1_path)
lib_lcmode1.gen_Teff_distribution.argtypes = [
    ctypes.POINTER(ctypes.c_float), ctypes.c_uint32,
    ctypes.c_float, ctypes.c_float, 
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.c_float, ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_float)
]
lib_lcmode1.gen_Teff_distribution.restype = None
lib_lcmode1.gen_lc.argtypes = [
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float), ctypes.c_float,
    ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), 
    ctypes.c_float, ctypes.c_float,
    ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_uint32), ctypes.c_uint32, ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_float), ctypes.c_uint32
]
lib_lcmode1.gen_lc.restype = None

