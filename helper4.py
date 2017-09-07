import numpy as np
import matplotlib.pyplot as plt

from setutil import PARAMS,Proton

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
    z1 = z+zpjump
    dz = np.mod(z1,2*zpjump)
    z1 = -zpjump+dz
    phix = factor*zpjump          # phase jump
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

def test0():
    t = np.arange(0.0, 6*np.pi, 2*np.pi/100.)
    dphi=np.radians(-40.)
    plt.subplot(121)
    plt.plot(t, f(t), 'b-') 
    plt.plot(t, f(t-dphi), 'g-')
    plt.plot(t, g(t,dphi), 'r-')
def test1():
    particle = Proton(tkin=100.)
    beta     = particle.beta
    lamb     = PARAMS['wellenlänge']
    gap      = PARAMS['spalt_laenge']
    phis     = [0,-25,-50,-75]
    for cnt,dphi in enumerate(phis):
        dphi = np.radians(dphi)
        zpjump  = (gap/2.+0.002)      # phase jump @ loc. of cavity ext. limit
        sig = zpjump/3.               # 3 sigma field strength @ cavity join
        mu1 = 0.                      # center of 1st cav. @ z=0
        mu2 = 2*zpjump                # center of 2nd cavity
        bl    = beta*lamb             # beta*lambda factor
        zr1  = mu2+zpjump             # right limit of intervall
        z  = np.arange(-zpjump,zr1,(zr1-zpjump)/100.)
        E  = [Ez(x,sig,mu1,mu2,dphi,bl,zpjump)[0] for x in z] # what the particle sees
        RF = [Ez(x,sig,mu1,mu2,dphi,bl,zpjump)[1] for x in z] # the time dependant modulation of the cavity field
        Ez0 =[A(x,sig,mu1,mu2) for x in z]                    # E(z,r=0) in cavities
        step = z[1]-z[0]
        Eint = 0.
        for i in range(1,len(z)):
            Eint += E[i]
        Eint = Eint*step/(z[len(z)-1]-z[0])
        Ezav = [Eint for x in z]      # average of field the particle sees
        
        ax = plt.subplot(2,2,cnt+1)
        ax.plot(z,E,   'b-',  label='Ez acc.(z)')
        ax.plot(z,RF,  'r-',  label=r'cos($phi$(z))')
        ax.plot(z,Ez0, 'g-',  label='Ez cav.(z)')
        ax.plot(z,Ezav,'k--', label='<Ez0> acc.')
        ax.set_title('sync.phase {:5.1f}[deg]'.format(np.degrees(dphi)))
        plt.legend(loc='lower right',fontsize='x-small')
    # ax.annotate('RF',           xy=(10,  0.74),   xytext=(14, 0.74),   arrowprops=dict(facecolor='red',   shrink=0.001))
    # ax.annotate('cav. profile', xy=(1.57,0.91),   xytext=(7.3, 0.91),  arrowprops=dict(facecolor='green', shrink=0.001))
    # ax.annotate('acc. field',   xy=(17.43, 0.56), xytext=(21.33, 0.56),arrowprops=dict(facecolor='blue',  shrink=0.001))
    # ax.annotate('Ez0 av.',      xy=(0.20, 0.32),  xytext=(0.20, 0.10), arrowprops=dict(facecolor='black', shrink=0.001))
    plt.show()

if __name__ == '__main__':
    # test0()
    test1()
