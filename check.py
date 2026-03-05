import multireflects
import numpy as np
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
import time
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

phi = np.linspace(0.0, 1.0, 500, dtype=np.float32)
lc_obj = multireflects.lc(phi, 1)
lc_obj.T1 = 25000.0; lc_obj.T2 = 3000.0
lc_obj.r1 = 0.3; lc_obj.r2 = 0.2
lc_obj.A1 = 0.9; lc_obj.A2 = 0.5
lc_obj.gdc1 = 1.0; lc_obj.gdc2 = 0.32
lc_obj.u1 = 0.1; lc_obj.u2 = 0.1
lc_obj.incl = 60.0; lc_obj.q = 1.0
lc_obj.P = 0.15; lc_obj.sma = 3.0
lc_obj.ecc = 0.0; lc_obj.w2 = 0.0
lc_obj.RV_gamma = 0.0

path = os.path.join(BASE_DIR, "TESS_Red_tcurve.dat")
S_curve_data = np.loadtxt(path, dtype=np.float32)
wave = S_curve_data[:,0]
S_curve = S_curve_data[:,1]

logT = np.linspace(np.log10(2000.0), 6.0, 100, dtype=np.float32)
T = 10.0**logT

wave_brd = wave[:, np.newaxis]
wave_brd = np.broadcast_to(wave_brd, (len(wave), len(logT)))
T_brd = T[:, np.newaxis]
T_brd = np.broadcast_to(T_brd, (len(T), len(wave)))
T_brd = np.transpose(T_brd, (1,0))
I_lambda = 2.0*const.h*const.c**2.0/(wave_brd*u.AA)**5.0 / (np.exp(const.h*const.c/const.k_B/(T_brd*u.K)/(wave_brd*u.AA)) - 1.0)
I_lambda = I_lambda.to(u.erg/u.s/(u.cm*u.cm)/u.AA).value
I_lambda = I_lambda.T

passband_obj = multireflects.passband(T, wave, I_lambda, S_curve, I_lambda)
passband_obj.u1 = np.array([0.0, 0.1, 0.0, 0.0], dtype=np.float32)
passband_obj.u2 = np.array([0.0, 0.1, 0.0, 0.0], dtype=np.float32)
passband_obj.A1 = 0.9
passband_obj.A2 = 0.5

passbandlist = [passband_obj,]

starttime = time.time()
for i in range(100):
    lc_obj.lc_mode0(passbandlist)
stoptime = time.time()
timespent = stoptime - starttime
print(f"{timespent:.3f} s")

lc = lc_obj.lc[0] / np.mean(lc_obj.lc[0])
plt.plot(phi, lc, linewidth=1)
path = os.path.join(BASE_DIR, "test.png")
plt.savefig(path, format="png", dpi=300)

