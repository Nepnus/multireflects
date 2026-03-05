import numpy as np

def get_shape(q, sma, b):
    RL= 0.49/(q**(2.0/3.0))*sma / (0.6/(q**(2.0/3.0)) + np.log(1.0+(1.0/q)**(1.0/3.0)))
    a = b; c = b; b2 = b
    tol = 1e-10
    phi_b = -sma**2.0/np.sqrt(b**2.0+sma**2.0) - q*sma**2.0/b
    
    while True:
        d_phi_a = -sma*2.0 / (sma - a)**2.0 + q*sma**2.0 / a**2.0 + 1.0 - (1.0+q)*a/sma
        phi_a = -sma**2.0/(sma-a) - q*sma**2.0/a + a - (1.0+q)/2.0/sma*a**2.0
        phi_ab = phi_a - phi_b
        a_new = a - phi_ab / d_phi_a
        if a_new >= b and a_new <= RL:
            if abs(a_new - a) < tol:
                a = a_new
                break
            a = a_new
        else:
            a = None
            break
        
    while True:
        d_phi_c = sma*2.0 / (sma + c)**2.0 + q*sma**2.0 / c**2.0 + 1.0 - (1.0+q)*c/sma
        phi_c = -sma**2.0/(sma+c) - q*sma**2.0/c - c - (1.0+q)/2.0/sma*c**2.0
        phi_cb = phi_c - phi_b
        c_new = c - phi_cb / d_phi_c
        if c_new >= b and c_new <= RL:
            if abs(c_new - c) < tol:
                c = c_new
                break
            c = c_new
        else:
            c = None
            break
        
    while True:
        d_phi_b2 = sma**2.0*b2/(b2**2.0+sma**2.0)**1.5 + q*sma**2.0/b2**2.0 - (1.0+q)*b2/sma
        phi_b2 = -sma**2.0/np.sqrt(sma**2.0+b2**2.0) - q*sma**2.0/b2 - (1.0+q)/2.0/sma*b2**2.0
        phi_b2b = phi_b2 - phi_b
        b2_new = b2 - phi_b2b / d_phi_b2
        if b2_new >= b and b2_new <= RL:
            if abs(b2_new - b2) < tol:
                b2 = b2_new
                break
            b2 = b2_new
        else:
            b2 = None
            break
        
    g_abratio = None; g_cbratio = None; g_b2bratio = None
    if a is not None and c is not None and b2 is not None:
        ga = sma**3.0 / (1.0+q) * (q/a**2.0 - 1.0/(sma-a)**2.0) + sma/(1.0+q) - a
        gc = sma**3.0 / (1.0+q) * (q/c**2.0 + 1.0/(sma+c)**2.0) - sma/(1.0+q) - c
        gb = sma**3.0 / (1.0+q) * (q/b**2.0 + b/(sma**2.0 + c**2.0)**1.5)
        gb2 = sma**3.0 / (1.0+q) * (b2/(b2**2.0 + sma**2.0)**1.5 + q/b2**2.0) - b2
        g_abratio = ga / gb
        g_cbratio = gc / gb
        g_b2bratio = gb2 / gb
    else:
        raise Exception("Roche lobe overflow.")
    
    return a, b2, c, g_abratio, g_b2bratio, g_cbratio

