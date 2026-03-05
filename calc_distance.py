import numpy as np

def keplerEquation_solve(M, ecc, tol=1e-10, max_iter=100):
    E_mid = None; E_low = None; E_high = None
    if ecc <= 0.65:
        sin_M = np.sin(M); sin_2M = np.sin(2.0*M); sin_3M = np.sin(3.0*M)
        E_mid = M + ecc*sin_M + 0.5*ecc**2.0*sin_2M + 0.125*ecc**3.0*(3.0*sin_3M - sin_M)
        for _ in range(max_iter):
            f = E_mid - ecc * np.sin(E_mid) - M
            f_prime = 1 - ecc * np.cos(E_mid)
            delta = f / f_prime
            E_mid -= delta
            if np.all(np.abs(delta) < tol):
                break
    else:
        E_low = np.zeros_like(M)
        E_high = 2.0*np.pi + E_low
        E_mid = (E_low + E_high) / 2
        for _ in range(max_iter):
            f_mid = E_mid - ecc*np.sin(E_mid) - M
            f_low = E_low - ecc*np.sin(E_low) - M
            left_mask = f_mid * f_low < 0
            E_high[left_mask] = E_mid[left_mask]
            E_low[~left_mask] = E_mid[~left_mask]
            E_mid = (E_low + E_high) / 2
            if np.all(E_high - E_low < tol):
                break
    E_mid[E_mid<0.0] += 2.0*np.pi
    E_mid[E_mid>=2.0*np.pi] -= 2.0*np.pi
    return E_mid

def calc_d(phi, ecc, w):
    v_init = np.pi/2.0 - w
    if v_init < 0.0:
        v_init += 2.0*np.pi
    E_init = 2.0*np.arctan2(np.sqrt(1 - ecc) * np.sin(v_init / 2), np.sqrt(1 + ecc) * np.cos(v_init / 2))
    M_init = E_init - ecc * np.sin(E_init)
    M = 2.0*np.pi*phi + M_init
    ind = M >= 2.0*np.pi
    M[ind] = M[ind] - 2.0*np.pi
    E = keplerEquation_solve(M, ecc)
    d = 1.0 - ecc*np.cos(E)                     # unit of d is 1.0/sma
    
    v = 2.0*np.arctan2(np.sqrt(1 + ecc) * np.sin(E / 2.0), np.sqrt(1 - ecc) * np.cos(E / 2.0))
    nv = v - (np.pi/2.0 - w)
    ind = nv >= 2.0*np.pi
    nv[ind] = nv[ind] - 2.0*np.pi
    ind = nv < 0.0
    nv[ind] = nv[ind] + 2.0*np.pi
    return d, nv

def calc_rv(nv, ecc, w, incl, q, sma, P):   # w rad, incl deg, q m2/m1, sma solar radius, P days
    P_s = P * 24.0 * 3600.0
    sma_km = sma * 695500.0
    incl_rad = incl * np.pi / 180.0
    
    v = nv + (np.pi/2.0 - w)
    ind = v >= 2.0*np.pi
    v[ind] = v[ind] - 2.0*np.pi
    ind = v < 0.0
    v[ind] = v[ind] + 2.0*np.pi
    
    RV2 = np.sin(incl_rad) * 2.0*np.pi*sma_km/P_s/(1.0+q) * (np.cos(v+w)+ecc*np.cos(w))/np.sqrt(1-ecc**2.0)
    
    w_1 = np.pi + w
    if w_1 >= 2.0*np.pi:
        w_1 -= 2.0*np.pi
    RV1 = np.sin(incl_rad) * 2.0*np.pi*sma_km/P_s/(1.0+1.0/q) * (np.cos(v+w_1)+ecc*np.cos(w_1))/np.sqrt(1-ecc**2.0)
    
    return RV1, RV2

def convertParams(fs, fc, star2flag=True):
    ecc = fs**2.0 + fc**2.0
    w = 0.0
    if ecc > 1.0e-15:
        w = np.arctan2(fs, fc)
        if w < 0.0:
            w += 2.0*np.pi
        if star2flag == False:
            w = np.pi + w
            if w >= 2.0*np.pi:
                w -= 2.0*np.pi
    return ecc, w

