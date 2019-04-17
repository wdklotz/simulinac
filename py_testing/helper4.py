import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0,'..')
from setutil import PARAMS,Proton,DEBUG

def NGauss(x,sigma,mu):    # Gauss Normalverteilung
    res = np.exp(-(((x-mu)/sigma)**2/2.))
    return res
def A(z,sigma,mu1,mu2):    # Amplitudenfunktion von 2 hintereinander liegenden Kavitaeten
    res = NGauss(z,sigma,mu1)+NGauss(z,sigma,mu2)
    return res
def Ez(z,sigma,mu1,mu2,dphi,bl,zpjmp):
    factor = 2.*np.pi/bl
    z0 = z                      # Ordinate
    Amp = A(z0,sigma,mu1,mu2)
    phix = factor*zpjmp         # phase jump
    phi0 = factor*z0+dphi
    phi1 = phi0-2.*phix
    rf0  = np.cos(phi0)
    rf1  = np.cos(phi1)
    res = NGauss(z0,sigma,mu1)*rf0 + NGauss(z0,sigma,mu2)*rf1 # superposition of 2 cavities
    rf = rf0 if z <= zpjmp else rf1
    return (res,rf)
def Intg(sigma,bl,dphi): # use integral fomula 3.896.2.pp.480 from I.S.Gradshteyn
    ex = np.exp(-2*(np.pi*sigma/bl)**2)
    res = np.sqrt(2*np.pi)*sigma*ex*np.cos(dphi)  # Note: proportional to cos(dphi)! simple!
    return res

def test0():
    particle = Proton(tkin=100.)
    beta     = particle.beta
    freq     = 800.*1.e6    # Hz
    lamb     = PARAMS['clight']/freq
    gap      = 0.044
    # phis     = [0, -30, -60., -90.]
    phis     = [0,-25, -50., -75.]
    for cnt,dphi in enumerate(phis):
        dphi = np.radians(dphi)
        # sigma  = gap/2./3.          # 3 sigma field strength @ cavity join
        sigma  = gap/2./2.2          # our data
        sigma  = gap/2./1.89         # wie in ez0.py
        bl   = beta*lamb             # beta*lambda factor
        mu1  = 0.                    # center of 1st cav. @ z=0
        mu2  = gap                   # center of 2nd cavity
        # zpjmp= (gap/2.+0.005)      # loc. z of phase jump (@ ext. limit of cavity)
        zpjmp= (gap/2.+0.00)         # loc. z of phase jump (@ ext. limit of cavity)
        zr   = mu2+zpjmp            # right limit of intervall
        zl   = -zpjmp               # left limit of intervall
        step = (zr-zl)/1000.
        z    = np.arange(zl,zr,step)
        E    = [Ez(x,sigma,mu1,mu2,dphi,bl,zpjmp)[0] for x in z] # what the particle sees
        RF   = [Ez(x,sigma,mu1,mu2,dphi,bl,zpjmp)[1] for x in z] # the time dependant modulation of the cavity field
        Ez0  = [A(x,sigma,mu1,mu2) for x in z]                   # E(z,r=0) peak field in cavities

        Vz = []; sum = 0
        for i in range(0,len(z)):
            sum += 20*E[i]*step
            Vz.append(sum)   # acc. voltage over intervall in [50KV] units to fit vert. axis
        
        ax = plt.subplot(2,2,cnt+1)
        ax.plot(z,E,  'b-',  label='Ez[MV/m]')
        ax.plot(z,RF, 'r-',  label=r'cos($phi$(z))')
        ax.plot(z,Ez0,'g-',  label='Ez-cav.[Mv/m]')
        ax.plot(z,Vz, 'k--', label='Vz[50KV]')
        ax.set_title(r'f={:5.1f}[MHz] gap={:3.0f}[mm] gap/2$\sigma$={:3.1f} sync.phase={:5.1f}[deg] Tkin={:3.0f}[MeV]'.
            format(freq*1.e-6,gap*1.e3,zpjmp/sigma,np.degrees(dphi),particle.tkin))
        ax.axhline(linestyle=':',color='m')
        ax.legend(loc='lower right',fontsize='x-small')
    plt.show()

if __name__ == '__main__':
    test0()
