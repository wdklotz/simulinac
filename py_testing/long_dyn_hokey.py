"""
Hokeystick simulation with symplectic integrator

Hamiltonian: von T.Wangler uses closure and decorates (returns) 4 function pointers.
Gap: calls Hamiltonian to preset the partial derivatives by closure. The returned 
    function pointers of the conditioned derivatives are uses by the
    integrator to calculate the mapping. Gap uses slices (nbslices) to model a single 
    gap-transition. Each slice presents a single interation step to map from 
    slice-input to slice-output.
Oder3Map: the 3.oder symplectic integrator from R.Ruth which makes a single interation
    step.

W.D.Klotz Mars 2019
"""
from math import pi,cos,sin,radians,degrees,sqrt
import scipy.constants as C
import matplotlib.pyplot as plt
import random as rn
    
def Hamiltonian(Ts,qE0T,phis,lamb):       
    """ 
        kanonische Koordinaten: p (aka w), x (aka phi)
        unabhaegige Variable: t (aka s)
        Hamiltonian H(p,x,t) = g(p) + V(x,t)
        g = A*p**2/2
        V = B*(sin(x)-x*cos(phis))
        f = -dV/dx = -B*(cos(x)-cos(phis))
        dgdp = dH/dp = A*p
    """
    # closure
    m0c2 = C.value('proton mass energy equivalent in MeV')  # MeV
    g    = 1+Ts/m0c2
    bgs  = sqrt(g**2-1.)
    A    = 2.*pi/bgs**3/lamb
    B    = qE0T/m0c2
    # potential
    def V(x,B=B,phis=phis):
        return B*(sin(x)-x*cos(phis))
    # hamiltonian
    def H(p,x,B=B,phis=phis):
        return A*p**2+V(x)
    # dH/dp
    def dgdp(p,A=A):
        return A*p
    # -dV/dx
    def f(x,t,B=B,phis=phis):
        return -B*(cos(x)-cos(phis))

    return dgdp,f,V,H

class Gap(object):
    def __init__(self,Ts,qE0T,phis,lamb,length=22.e-3):
        self.Ts       = Ts
        self.Tsf      = Ts     # final kin. energy
        self.qE0T     = qE0T
        self.phis     = phis   # rad
        self.lamb     = lamb
        self.length   = length
        self.nbslices = 10     # nbof slices
        dl            = length/self.nbslices
        m0c2          = C.value('proton mass energy equivalent in MeV')  # MeV
        g             = 1.+Ts/m0c2   # Lorentz gamma
        bg            = sqrt(g**2-1.)
        beta          = bg/g         # Lorentz beta
        # phases,energies & maps of slices
        phases    = []
        dtks      = []
        self.maps = []
        for i in range(self.nbslices+1):
            phase = self.phis+2*pi/(beta*lamb)*(i*dl-self.length/2.)
            phases.append(phase)
            dtk = qE0T*dl*cos(phase)
            dtks.append(dtk)
        for i in range(self.nbslices):
            dgdp,f,pot,ham = Hamiltonian(self.Tsf,qE0T,phases[i],lamb)
            self.maps.append(Order3Map(+dl,dgdp,f))
            self.Tsf += dtks[i]
        self.Tsf -= dtks[i]     # one too much
    def map(self,v):
        cnt=0
        while cnt < self.nbslices:
            v = self.maps[cnt].map(v)
            cnt += 1
        return v
        # init again with different kin. energy
    def __call__(self,Ts):
        self.__init__(Ts,self.qE0T,self.phis,self.lamb)
        
class Order3Map(object):
    """
    Ruth IEEE Transactions on Nuclear Science, Vol. NS-30, No. 4, August 1983
    """
    def __init__(self,stepsize,dgdp,f):
        self.h    = stepsize
        self.c1   = 7./24.
        self.d1   = 2./3.
        self.c2   = 3./4.
        self.d2   = -2./3.
        self.c3   = -1./24.
        self.d3   = 1.
        self.gp   = dgdp     # dg/dp
        self.f    = f        # -dV/dx
    def map(self,v0):        # phase vector: v=(p,x,t)
        v1 = self.step1(v0)
        v2 = self.step2(v1)
        v3 = self.step3(v2)
        return v3
    def step1(self,v):
        p0,x0,t0 = v
        p1 = p0 + self.c1 * self.h * self.f(x0,t0)
        x1 = x0 + self.d1 * self.h * self.gp(p1)
        return (p1,x1,t0)
    def step2(self,v):
        p1,x1,t0 = v
        p2 = p1 + self.c2 * self.h * self.f(x1,t0+self.d1*self.h)
        x2 = x1 + self.d2 * self.h * self.gp(p2)
        return (p2,x2,t0)
    def step3(self,v):
        p2,x2,t0 = v
        p = p2 + self.c3 * self.h * self.f(x2,t0+(self.d1+self.d2)*self.h)
        x = x2 + self.d3 * self.h * self.gp(p)
        t = t0 + self.h
        return (p,x,t)
    
#todo: besser loop particles in loop gaps
def test1(Ts,qE0T,phis,lamb):
    print("----------------------Test 1--------------")
    # prep plot
    m0c2   = C.value('proton mass energy equivalent in MeV')  # MeV
    phimin = radians(-60.)
    phimax = -phimin
    # dphi = (phimax-phimin)/20.
    wmin = -0.05/m0c2
    wmax = +0.05/m0c2
    # dw = (wmax-wmin)/20.
    pou    = []
    pin    = []
    nbprtcls  = 5000
    nbgaps    = 40
    # loop particles
    for j in range(nbprtcls):
        # x0 = phi0 = rn.uniform(phimin,phimax)
        # p0 = w0   = rn.uniform(wmin,wmax)
        x0 = phi0 = rn.gauss(phis,2*phis)
        p0 = w0   = rn.gauss(0,2*wmax)
        t0 = 0
        v = (p0,x0,t0)
        pin.append(v)
        gap  = Gap(Ts,qE0T,phis,lamb)
        # loop gaps
        for i in range(nbgaps):
            v = gap.map(v)
            # increase particle impulse for next gap
            gap(gap.Tsf)
        pou.append(v)
    print("Ts {}".format(gap.Tsf))
 
    # plot results
    win   = [x[0]*m0c2     for x in pin]
    phin  = [degrees(x[1]) for x in pin]
    wf    = [x[0]*m0c2     for x in pou]
    phif  = [degrees(x[1]) for x in pou]
 
    ax1   = plt.subplot(121)
    plt.xlim([-100,200])
    plt.ylim([-0.4,1.2])
    ax1.autoscale(enable=False)
    ax1.scatter(phin,win,s=1,c='r')

    ax2   = plt.subplot(122,sharex=ax1,sharey=ax1)
    ax2.autoscale(enable=False)
    ax2.scatter(phif,wf,s=1,c='r')
    plt.show()
    
def test0(Ts,qE0T,phis,lamb):
    print("----------------------Test 0--------------")
    m0c2 = C.value('proton mass energy equivalent in MeV')  # MeV
    dgdp,dfdx,pot,ham  = Hamiltonian(Ts,qE0T,phis,lamb)
    dt = 1.e-2
    map3p = Order3Map(+dt,dgdp,dfdx)
    map3m = Order3Map(-dt,dgdp,dfdx)

    ax = plt.subplot()
    phimin = 1.5*phis
    phimax = -phis
    dphi = (phimax-phimin)/20.
    phi0 = phimin

    wmin = -0.05/m0c2
    wmax = +0.05/m0c2
    dw = (wmax-wmin)/20.
    w0 = wmin
    # while phi0 <= phimax:
    for i in range(50):
        # phi0 = rn.gauss(degrees(phis),1.*degrees(abs(phis)))
        phi0 = rn.uniform(phimin,phimax)
        w0   = rn.uniform(wmin,wmax)

        t0 = 0.
        p0 = w0
        x0 = phi0
        vp = (p0,x0,t0)
        vm = (p0,x0,t0)
        pou = [vp]
        for step in range(300):
            vp = map3p.map(vp)
            vm = map3m.map(vm)
            pou.append(vp)
            pou.append(vm)
        # phi0 += dphi
        w   = [x[0]*m0c2     for x in pou]
        phi = [degrees(x[1]) for x in pou]
        plt.xlim([-60,60])
        plt.ylim([-0.1,0.1])
        plt.autoscale(enable=False,axis='both')
        ax.scatter(phi,w,s=1,c='r')
    plt.show()

if __name__ == '__main__':
    clight        = C.c          # m/s
    freq          = 408.e6       # Hz
    lamb          = clight/freq  # m
    qE0T          = 1.           # MeV/m
    phis          = radians(-20.)
    Ts            = 5.
    
    test0(Ts,qE0T,phis,lamb)    
    test1(Ts,qE0T,phis,lamb)