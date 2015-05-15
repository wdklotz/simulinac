#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from pylab import plot,show,legend
from math import cos,sin,radians,degrees
'''try Ruth's symplectic integrator'''
def display(functions):
    for function in functions:
        t = [x[0] for x in function]
        p = [x[1] for x in function]
        q = [x[2] for x in function]
        vp= [x[3] for x in function]
        # plot(t,p,label='p')
        # plot(t,q,label='q')
        # legend(loc='lower right',fontsize='x-small')
        plot(q,p,label='w(phi)',color='blue')
        plot(q,vp,label='vp(phi)',color='red')
    show()
class Order3Map(object):
    def __init__(self,stepsize,dHdP,dHdX):
        self.h = stepsize
        self.c1 = 7./24.
        self.d1 = 2./3.
        self.c2 = 3./4.
        self.d2 = -2./3.
        self.c3 = -1./24.
        self.d3 = 1.
        self.g  = dHdP     # dH/dp
        self.f  = dHdX     # dH/dx
    def map(self,v0):
        v1 = self.step1(v0)
        v2 = self.step2(v1)
        v3 = self.step3(v2)
        return v3
    def step1(self,v):
        p0 = v[0]
        x0 = v[1]
        p1 = p0 - self.c1 * self.h * self.f(x0)
        x1 = x0 + self.d1 * self.h * self.g(p1)
        return (p1,x1)
    def step2(self,v):
        p1 = v[0]
        x1 = v[1]
        p2 = p1 - self.c2 * self.h * self.f(x1)
        x2 = x1 + self.d2 * self.h * self.g(p2)
        return (p2,x2)
    def step3(self,v):
        p2 = v[0]
        x2 = v[1]
        p3 = p2 - self.c3 * self.h * self.f(x2)
        x3 = x2 + self.d3 * self.h * self.g(p3)
        return (p3,x3)
def dHdP(p):   
    '''canonical equation dH/dp = dx/dt == p; H = p**2/2+V(x)'''
    return p
def dHdX(x):
    '''canonical equation dH/dx = - dp/dt = - dV/dx for Pendel Hamiltonian V(x)=cos(x)'''
    return +sin(x)
def dTWdX(cphis):
    '''canonical equation dH/dx = - dV/dx; TWrangler's' Hamiltonian V(x)=sin(x)-x*cos(x0))'''
    return lambda x: +(cos(x)+cphis)   # closure!    
def Vp(phi,phis):
    '''pseudo potential in TW's Hamiltonian V(x)=sinc(x)-x*cos(x0)'''
    return +(sin(phi)-phi*cos(phis))
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
def test0(): 
    print(
    '''Shows symplectic integration of:
    Hamiltonian for TWrangler's gap (similar to Lapostolle CERN 87-09)
    H = p**2/2 +V(x) = p**2/2 + (sin(x) - x * cos(phis))
    p stands for (scaled) w and V(x) for V(phase phi) w/o dimensional scaling and
    overlays a plot of V(x)''')
    h = 1e-2              # integrator's step size'
    phis = radians(-90)   # parameter cos(phi0)
    cphis = cos(phis)
    map = Order3Map(h,dHdP,dTWdX(cphis))
    q0   = radians(-90)   # anfangswert q0
    # pmax=+1.4134
    pmax=+3.0005          # max p
    anz = 30              # anz trajectories
    dp = pmax/anz
    pmin=-10*dp
    pvalues = [pmin+i*dp for i in range(anz+1)]  # list of anfangswerte p0
    functions=[]
    for pv in pvalues:  # loop over several p
        v = (pv,q0)
        function=[]
        for i in range(int(10./h)):   # loop over independent variable (time z.B.)
            pot = Vp(v[1],phis)
            function.append((i,v[0],degrees(v[1]),pot))
            v = map.map(v)            # [p(t),q(t)] <---map---- [p(0),q(0)]
        functions.append(function)
    display(functions)
    print(
    '''
    =============================================================
    I still don't understand the outcome of this excercise???!!!!
    =============================================================''')
def test1(): 
    print(
    '''Shows symplectic integration of:
    Hamiltonian for TWrangler's gap (similar to Lapostolle CERN 87-09)
    H = p**2/2 +V(x) = p**2/2 + (sin(x) - x * cos(phis))
    p stands for (scaled) w and V(x) for V(phase phi) w/o dimensional scaling and
    overlays a plot of V(x)''')
    h = 1e-3                # integrator's step size'
    # phis = radians(-90)   # parameter cos(phi0)
    # cphis = cos(phis)
    map = Order3Map(h,dHdP,dHdX)
    q0   = radians(0)   # anfangswert q0
    pmax=2
    pmin=1e-3
    anz = 50            # anz trajectories
    dp = (pmax-pmin)/anz
    pvalues = [pmin+i*dp for i in range(0,anz+1,3)]  # list of anfangswerte p0
    functions=[]
    for pv in pvalues:  # loop over several p
        v = (pv,q0)
        function=[]
        for i in range(int(28./h)):   # loop over independent variable (time z.B.)
            pot = -cos(v[1])
            function.append((i,v[0],degrees(v[1]),pot))
            v = map.map(v)            # [p(t),q(t)] <---map---- [p(0),q(0)]
        functions.append(function)
    display(functions)
    print(
    '''
    =============================================================
    I still don't understand the outcome of this excercise???!!!!
    =============================================================''')
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    # test0()
    test1()

    
    
    
