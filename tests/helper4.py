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
def Ez(z,sig,mu1,mu2,dphi,bl,zpjump):
    factor = 2.*np.pi/bl
    z0 = z                      # Ordinate
    Amp = A(z0,sig,mu1,mu2)
    phix = factor*zpjump        # phase jump
    phi0 = factor*z0+dphi
    phi1 = phi0-2.*phix
    rf0  = np.cos(phi0)
    rf1  = np.cos(phi1)
    if z <= zpjump:
        rf = rf0          # in 1st cavity
    else:
       rf = rf1           # in 2nd cavity
    res = NGauss(z0,sig,mu1)*rf0 + NGauss(z0,sig,mu2)*rf1 # superposition of 2 cavities
    return (res,rf)
def Intg(sig,bl,dphi): # use integral fomula 3.896.2.pp.480 from I.S.Gradshteyn
    ex = np.exp(-2*(np.pi*sig/bl)**2)
    res = np.sqrt(2*np.pi)*sig*ex*np.cos(dphi)  # Note: proportional to cos(dphi)! simple!
    return res

def test0():
    particle = Proton(tkin=100.)
    beta     = particle.beta
    lamb     = PARAMS['wellenlänge']
    gap      = PARAMS['spalt_laenge']
    phis     = [0,-25,-50., -75.]
    for cnt,dphi in enumerate(phis):
        dphi = np.radians(dphi)
        zpjump  = (gap/2.+0.005)      # loc. of phase jump @ ext. limit of cavity 
        sig  = zpjump/3.              # 3 sigma field strength @ cavity join
        mu1  = 0.                     # center of 1st cav. @ z=0
        mu2  = 2*zpjump               # center of 2nd cavity
        bl   = beta*lamb              # beta*lambda factor
        zr   = mu2+zpjump             # right limit of intervall
        zl   = -zpjump                # left limit of intervall
        step = (zr-zl)/1000.
        z    = np.arange(zl,zr,step)
        E    = [Ez(x,sig,mu1,mu2,dphi,bl,zpjump)[0] for x in z] # what the particle sees
        RF   = [Ez(x,sig,mu1,mu2,dphi,bl,zpjump)[1] for x in z] # the time dependant modulation of the cavity field
        Ez0  = [A(x,sig,mu1,mu2) for x in z]                    # E(z,r=0) in cavities

        Ez_av_int = 2*Intg(sig,bl,dphi)/(zr-zl)   # average using integral formula
        Ezav = [Ez_av_int for x in z]             # average of field the particle sees

        # Ezsum = 0.                              # numerical integration
        # for i in range(0,len(z)):
        #     Ezsum += E[i]
        # Ez_av_num = Ezsum*step/(zr-zl)
        # DEBUG('<Ez>num {}  <Ez>int {}'.format(Ez_av_num,Ez_av_int))
        
        ax = plt.subplot(2,2,cnt+1)
        ax.plot(z,E,   'b-',  label='Ez acc.(z)')
        ax.plot(z,RF,  'r-',  label=r'cos($phi$(z))')
        ax.plot(z,Ez0, 'g-',  label='Ez cav.(z)')
        ax.plot(z,Ezav,'k--', label='<Ez0> acc.')
        ax.set_title('sync.phase{:5.1f}[deg], T(p+){:5.1f}[MeV]'.format(np.degrees(dphi),particle.tkin))
        plt.legend(loc='lower right',fontsize='x-small')
    plt.show()

if __name__ == '__main__':
    test0()
