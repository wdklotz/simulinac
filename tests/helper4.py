import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0,'..')
from setutil import PARAMS,Proton,DEBUG

def NGauss(x,sig,mu):    # Gauss Normalverteilung
    res = np.exp(-(((x-mu)/sig)**2/2.))
    return res
def A(z,sig,mu1,mu2):    # Amplitudenfunktion von 2 hintereinander liegenden Kavitaeten
    res = NGauss(z,sig,mu1)+NGauss(z,sig,mu2)
    return res
def Ez(z,sig,mu1,mu2,dphi,bl,zpjmp):
    factor = 2.*np.pi/bl
    z0 = z                      # Ordinate
    Amp = A(z0,sig,mu1,mu2)
    phix = factor*zpjmp         # phase jump
    phi0 = factor*z0+dphi
    phi1 = phi0-2.*phix
    rf0  = np.cos(phi0)
    rf1  = np.cos(phi1)
    res = NGauss(z0,sig,mu1)*rf0 + NGauss(z0,sig,mu2)*rf1 # superposition of 2 cavities
    rf = rf0 if z <= zpjmp else rf1
    return (res,rf)
def Intg(sig,bl,dphi): # use integral fomula 3.896.2.pp.480 from I.S.Gradshteyn
    ex = np.exp(-2*(np.pi*sig/bl)**2)
    res = np.sqrt(2*np.pi)*sig*ex*np.cos(dphi)  # Note: proportional to cos(dphi)! simple!
    return res

def test0():
    particle = Proton(tkin=100.)
    beta     = particle.beta
    freq     = 800.*1.e6    # Hz
    lamb     = PARAMS['lichtgeschwindigkeit']/freq
    gap      = PARAMS['spalt_laenge']
    gap      = 0.044
    # phis     = [0, -30, -60., -90.]
    phis     = [0,-25, -50., -75.]
    for cnt,dphi in enumerate(phis):
        dphi = np.radians(dphi)
        zpjmp= (gap/2.+0.005)        # loc. z of phase jump (@ ext. limit of cavity) 
        # sig  = zpjmp/3.              # 3 sigma field strength @ cavity join
        sig  = zpjmp/2.2             # our data
        # sig  = zpjmp/1.              # test
        bl   = beta*lamb             # beta*lambda factor
        mu1  = 0.                    # center of 1st cav. @ z=0
        mu2  = 2*zpjmp               # center of 2nd cavity
        zr   = mu2+zpjmp             # right limit of intervall
        zl   = -zpjmp                # left limit of intervall
        step = (zr-zl)/1000.
        z    = np.arange(zl,zr,step)
        E    = [Ez(x,sig,mu1,mu2,dphi,bl,zpjmp)[0] for x in z] # what the particle sees
        RF   = [Ez(x,sig,mu1,mu2,dphi,bl,zpjmp)[1] for x in z] # the time dependant modulation of the cavity field
        Ez0  = [A(x,sig,mu1,mu2) for x in z]                   # E(z,r=0) peak field in cavities

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
            format(freq*1.e-6,gap*1.e3,zpjmp/sig,np.degrees(dphi),particle.tkin))
        ax.axhline(linestyle=':',color='m')
        ax.legend(loc='lower right',fontsize='x-small')
    plt.show()

if __name__ == '__main__':
    test0()
