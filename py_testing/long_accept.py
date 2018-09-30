import math as M
import scipy.constants as C

## Long. Emittance
def waccept(E0T,T,freq,phis):
    """
    Helper to calculate longitudinal phase space ellipse parameters nach T.Wangler (6.47-48) pp.185
        (w/w0)**2 + (Dphi/Dphi0)**2 = 1
        emitw = w0(x)Dphi0 = ellipse_area/pi
    IN
        E0T:   the gap field times ttf: E0 * TTF [Mv/m]
        T:     the kinetic energy of partice (here proton) [MeV]
        freq:  the rf frequency [Hz]
        phis:  the synchronous phase [rad]
    """

    pi             = C.pi
    lamb           = C.c/freq           # [m]
    m0c2,unit,prec = C.physical_constants['proton mass energy equivalent in MeV']
    gamma          = 1. + T/m0c2
    beta           = M.sqrt(1.-1./(gamma*gamma))
    gb             = gamma*beta

    # large amplitude oscillations (T.Wangler pp. 175)
    factor_phis = phis*M.cos(phis)-M.sin(phis)
    wmax  = M.sqrt(2.*E0T*gb**3*lamb/(pi*m0c2)*factor_phis)     # T.Wangler (6.28)
    DWmax = wmax*m0c2       # [MeV]
    phi_1 = -phis           # [rad]
    phi_2 = 2.*phis         # [rad] Naehrung T.Wangler pp.178
    psi   = 3.*M.fabs(phis) # [rad]
    # Dp2pmx  = gamma/(gamma+1)*DWmax/T*100.    # [%]
    Dp2pmx = gamma/(gamma*gamma-1)*wmax*100. # [%]

    res = (T,E0T,freq,wmax,DWmax,M.degrees(psi),Dp2pmx)
    return res
    
if __name__ == '__main__':
    E0T = 1
    Trange = [70.-i*5. for i in range(0,14)]
    freq = 800.e6
    phis = M.radians(-20.)
    
    for T in Trange:
        Tr, E0Tr, freqr, wmaxr, DWmaxr, psi, Dp2pmx = waccept(E0T,T,freq,phis)
        print('T {:.2f} [Mev] E0T {:.2f} [Mv/m] wx {:.2e} [rad] DW {:.2f} [MeV] Dp/p {:.2f} [%]'.format(Tr,E0Tr,wmaxr,DWmaxr,Dp2pmx))
    
