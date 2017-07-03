#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# -*- coding: utf-8 -*-
"""
Copyright 2016 Wolf-Dieter Klotz <wdklotz@gmail.com>
"""
from math import sqrt, sin, cos, sinh, cosh, log, acos, degrees
import matplotlib.pyplot as plt

def betagamma(tkin):    # 
    """
    Out: beta*gamma for proton
    In:
        tkin [Mev] kinetic energy
    """
    tkin_zu_e0 = tkin/928.2796
    res = 2.*tkin_zu_e0 + tkin_zu_e0**2
    res = sqrt(res)
    return res

def kq(grad,tkin):  
    """
    Out: kq [1/m^2] quad for proton 
    In:
        grad [T/m] gradient
        tkin [Mev] kin energy
    """
    res = 0.31952*grad/betagamma(tkin)
    return res

def bg(kq,tkin):    
    """
    Out: bg [T/m] gradient dB/ds for protons
    In:
        kq [1/m^2] quad
        tkin [Mev] kinetic energy
    """
    res = kq/0.31952*betagamma(tkin)
    return res

##-------------------------------------Test2--
def test0(params):
    """ 
    Wiedemann's matices for FODO with acceleration Vol II, chapter 6.4, pp 222
    """
    import elements
    elements.MDIM=2
    print('------------------------------Test0--')
    print('FODO Lattice and Acceleration [Wiedeman II,pp 220]')
    tkin = params['tkin']
    lq = params['lq']
    L = params['L']
    ldrift = L-lq
    lquad = lq/2.
    grad_f = params['grad_f']
    grad_d = params['grad_d']
    kfoc = kq(grad_f,tkin)
    # ssize = 0.04
    betgam = betagamma(tkin)
    gamma = 1.+ tkin/928.2796
    beta = betgam/gamma
    cp0 = 928.2796*betgam
    alpha = params['acc']   # acceleration alpha [MeV/m] Wiedeman II,(622.1),pp 222
    eta0 = alpha/beta/cp0
    print('Ldrift [m]  {:4.4f}, Lquad [m] {:4.4f}, Lhalf-cell [m] {:4.4f}'.format(ldrift,lq,L))
    print('grad_f [T/m] {:4.4f}, grad_d [T/m] {:4.4f}, kfoc [1/m^2] {:4.4f}'.format(grad_f,grad_d,kfoc))
    print('beta*gamma  {:4.4f}, gamma {:4.4f}, beta {:4.4f}'.format(betgam,gamma,beta))
    print('tk-in [Mev] {:4.4f}, cp0 [Mev] {:4.4f}'.format(tkin,cp0))
    print('alpha [MeV/m] {:4.4f} (Acceleration), eta0 [1/m] {:4.4f}'.format(alpha,eta0))

    # Mdrift Wiedeman II,(6.233)
    sigma4 = 1./(1.+eta0*ldrift)
    drift = elements.I(label='Drift')
    drift.matrix[0,1]=-log(sigma4)/eta0
    drift.matrix[1,1]=sigma4
    print(drift.string())
    
    def M(m,C,S):
        m.matrix[0,0]= sigma*C(dksi)+sigma/8.*(3./ksi0+1./2.)*S(dksi)
        m.matrix[0,1]= sigma/sqrt(k0)*S(dksi)+sigma/(8.*sqrt(k0))*dksi/(ksi0*ksil)*C(dksi)
        m.matrix[1,0]=-sigma**3*sqrt(k0)*S(dksi)+3.*(sigma**3)/8.*dksi/(ksi0*ksil)*sqrt(k0)*C(dksi)
        m.matrix[1,1]= sigma**3*C(dksi)-(sigma**3)/8.*(1./ksi0+3./2.)*S(dksi)
        return
        
    # Mfoc Wiedemann II,(6.234)
    sigma4=1./(1.+eta0*lquad)
    sigma = sqrt(sqrt(sigma4))
    k0 = kq(grad_f,tkin)
    ksi0 = 2.*sqrt(k0)/eta0
    ksil = 2.*sqrt(k0*(1.+eta0*lquad))/eta0
    dksi=ksil-ksi0
    # print('sigma4, sigma, ksi0, ksil,dksi {:4.4f} {:4.4f} {:4.4f} {:4.4f} {:4.4f}'.format(sigma4,sigma,ksi0,ksil,dksi))    
    print('approx(2*sqrt(k0)/eta0 >> 1) valid? {:4.4f}'.format(2.*sqrt(k0)/eta0))   
    print('approx(eta0*l/4        << 1) valid? {:4.4f}'.format(eta0*lquad/4.))
    mfoc = elements.I(label='Mfoc')
    M(mfoc,cos,sin)
    print(mfoc.string())
    
    # Mdefoc Wiedemann II,(6.237)
    k0 = -kq(grad_d,tkin)
    ksi0 = 2.*sqrt(k0)/eta0
    ksil = 2.*sqrt(k0*(1.+eta0*lquad))/eta0
    dksi=ksil-ksi0
    print('approx(2*sqrt(k0)/eta0 >> 1) valid? {:4.4f}'.format(2.*sqrt(k0)/eta0))
    print('approx(eta0*l/4        << 1) valid? {:4.4f}'.format(eta0*lquad/4.))
    mdfoc = elements.I(label='Mdfoc')
    M(mdfoc,cosh,sinh)
    print(mdfoc.string())
    
    # Mfodo
    fodo_cell = elements.I(label='I')
    fodo_cell = mfoc*fodo_cell
    fodo_cell = drift*fodo_cell
    fodo_cell = mdfoc*fodo_cell
    fodo_cell = mdfoc*fodo_cell
    fodo_cell = drift*fodo_cell
    fodo_cell = mfoc*fodo_cell
    fodo_cell.label = 'FODDOF'
    print(fodo_cell.string())
    
    trace = fodo_cell.matrix[0,0]+fodo_cell.matrix[1,1]
    stable = False
    if trace < 2.: stable = True
    print('1/2*Tr = {:4.4f} ==> stable = {}'.format(0.5*trace,stable))
    
    if stable:
        twbeta = fodo_cell.matrix[0,1]**2/(1.-fodo_cell.matrix[0,0]**2)
        twbeta = sqrt(twbeta)
        twgamma = 1./twbeta
        twalpha = 0.
        bphase = acos(fodo_cell.matrix[0,0])
        print('Twiss functions: alpha {:4.4f}, beta {:4.4f}, gamma {:4.4f}\nphase_advance [deg] {:4.4f}'
            .format(twalpha,twbeta,twgamma,degrees(bphase)))
#-------------------------------------------main---
if __name__ == '__main__':
    # wille test
    tk = 5.
    k = 1.2
    bgrad = bg(k,tk)
    wille=dict(
             L      = 3.    # 1/2 cell length
            ,lq     = 0.4   # full quad length
            ,tkin   = tk
            ,grad_f = +bgrad
            ,grad_d = -bgrad
            ,acc    = 1.e-1
            )
    test0(wille)