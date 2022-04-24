#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
"""
Copyright 2015 Wolf-Dieter Klotz <wdklotz@gmail.com>
This file is part of the SIMULINAC code

    SIMULINAC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    SIMULINAC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SIMULINAC.  If not, see <http://www.gnu.org/licenses/>.
"""
import matplotlib.pyplot as plt
from math import cos,sin,radians,degrees

'''try Ruth's symplectic integrator'''
def display(functions,text=''):
    for function in functions:
        t = [x[0] for x in function]
        p = [x[1] for x in function]
        q = [x[2] for x in function]
        vp= [x[3] for x in function]
        plt.figure(text)
        plt.plot(t,p,label='p')
        plt.plot(t,q,label='q')
        plt.legend(loc='lower right',fontsize='x-small')
        # plot(q,p,label='w(phi)',color='blue')
        # plot(q,vp,label='vp(phi)',color='red')
    plt.show()
class Order2Map(object):  # from Ruth IEEE Transactions on Nuclear Science, Vol. NS-30, No. 4, August 1983
    def __init__(self,stepsize,dHdP,dHdX):
        self.h = stepsize
        self.g  = dHdP     # dH/dp
        self.f  = dHdX     # dH/dx
    def map(self,v0):       # coordinate vector: v==(x,y,z)==(t,p,q)
        v1 = self.step1(v0)
        v2 = self.step2(v1)
        return v2
    def step1(self,v):
        x0 = v[0]
        y0 = v[1]
        z0 = v[2]
        z1 = z0
        y1 = y0 + 0.5 * self.h * z1
        return (x0,y1,z1)
    def step2(self,v):
        x1 = v[0]
        y1 = v[1]
        z1 = v[2]
        z2 = z1 + self.h * self.f(y1)
        y2 = y1 + 0.5 * self.h * z2
        x2 = x1 + self.h
        return (x2,y2,z2)
class Order3Map(object):  # from Ruth IEEE Transactions on Nuclear Science, Vol. NS-30, No. 4, August 1983
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
    def map(self,v0):       # coordinate vector: v==(x,y,z)==(t,p,q)
        v1 = self.step1(v0)
        v2 = self.step2(v1)
        v3 = self.step3(v2)
        return v3
    def step1(self,v):
        x0 = v[0]
        y0 = v[1]
        z0 = v[2]
        z1 = z0 + self.c1 * self.h * self.f(y0)
        y1 = y0 + self.d1 * self.h * self.g(z1)
        return (x0,y1,z1)
    def step2(self,v):
        x1 = v[0]
        y1 = v[1]
        z1 = v[2]
        z2 = z1 + self.c2 * self.h * self.f(y1)
        y2 = y1 + self.d2 * self.h * self.g(z2)
        return (x1,y2,z2)
    def step3(self,v):
        x2 = v[0]
        y2 = v[1]
        z2 = v[2]
        z3 = z2 + self.c3 * self.h * self.f(y2)
        y3 = y2 + self.d3 * self.h * self.g(z3)
        x3 = x2 + self.h
        return (x3,y3,z3)
class RKMap11_2(object):      # Runke-Kutta nach G.Jordan-Englen/F.Reutter (BI Hochschultaschenbücher Band 106 pp.240)
    def __init__(self,h,f,g):
        self.h = h
        self.f = f
        self.g = g
    def map(self,vi):       # coordinate vector: v==(x,y,z)==(t,q,p)
        xi = vi[0]
        yi = vi[1]
        zi = vi[2]
        h = self.h
        f = self.f
        g = self.g
        k1 = h * f(xi,yi,zi)
        l1 = h * g(xi,yi,zi)
        k2 = h * f(xi + 0.5*h, yi + 0.5*k1, zi + 0.5*l1)
        l2 = h * g(xi + 0.5*h, yi + 0.5*k1, zi + 0.5*l1)
        k3 = h * f(xi + 0.5*h, yi + 0.5*k2, zi + 0.5*l2)
        l3 = h * g(xi + 0.5*h, yi + 0.5*k2, zi + 0.5*l2)
        k4 = h * f(xi + h, yi + k3, zi +l3)
        l4 = h * g(xi + h, yi + k3, zi +l3)
        k = 1./6.*(k1 + 2.*k2 + 2.*k3 +k4)
        l = 1./6.*(l1 + 2.*l2 + 2.*l3 +l4)
        x = xi + h
        y = yi + k
        z = zi + l
        return (x,y,z)
class RKMap11_3(object):      # Runke-Kutta nach G.Jordan-Englen/F.Reutter (BI Hochschultaschenbücher Band 106 pp.240)
    def __init__(self,h,g):
        self.h = h
        self.g = g
    def map(self,vi):       # coordinate vector: v==(x,y,z)==(t,q,p)
        h = self.h
        g = self.g
        xi = vi[0]
        yi = vi[1]
        zi = vi[2]           # y'  = z
        y2p= g(xi,yi,zi)     # y'' = z' = g(x,y,z)
        l1 = h * y2p
        l2 = h * g(xi + 0.5*h, yi + 0.5*h*zi                , zi + 0.5*l1)
        l3 = h * g(xi + 0.5*h, yi + 0.5*h*zi + 0.25*h*h*y2p, zi + 0.5*l2)
        l4 = h * g(xi + h    , yi + h*zi + 0.5*h*l2        , zi + l3)
        x = xi + h
        y = yi + h*zi + h/6.*(l1 + l2 + l3)
        z = zi + 1./6.*(l1 + 2.*(l2 + l3) + l4)
        return (x,y,z)
def pendel_f(x,y,z):
    return z
def pendel_g(a=90.):
    c = cos(radians(a))
    def g(x,y,z):
        return -sin(y)-c
    return g
def dHdP(p):              #  dH/dp
    '''canonical equation dH/dp = dx/dt = p; H = p**2/2+V(x)'''
    return +p
def dHdX(x):              #  dH/dx
    '''canonical equation dH/dx = - dp/dt = - dV/dx for Pendel Hamiltonian V(x)=cos(x)'''
    return -sin(x)
def dTWdX(cphis):         #  dH/dx
    '''canonical equation dH/dx = - dV/dx; TWrangler's' Hamiltonian V(x)=sin(x)-x*cos(x0))'''
    return lambda x: -(cos(x)+cphis)   # closure!
def Vp(phi,phis):         #  V(x)
    '''pseudo potential in TW's Hamiltonian V(x)=sin(x)-x*cos(x0)'''
    return +(sin(phi)-phi*cos(phis))
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
def test0():              # use T.Wrangler
    print(
    '''Shows symplectic integration of:
    Hamiltonian for TWrangler's gap (similar to Lapostolle CERN 87-09)
    H = p**2/2 +V(q) = p**2/2 + (sin(q) - q * cos(phis))
    p stands for (scaled) w and V(q) for V(phase phi) w/o dimensional scaling and
    overlays a plot of V(q)''')
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
        v = (0.,pv,q0)
        function=[]
        for i in range(int(10./h)):   # loop over independent variable (time z.B.)
            pot = Vp(v[2],phis)
            function.append((i,v[1],degrees(v[2]),pot))
            v = map.map(v)            # [p(t),q(t)] <---map---- [p(0),q(0)]
        functions.append(function)
    display(functions,text='TW Order3Map')
    print(
    '''
    =============================================================
    I still don't understand the outcome of this excercise???!!!!
    =============================================================''')
def test1(order=3):              # use pendulum
    print(
    '''Shows symplectic integration of:
    Hamiltonian for pendulum
    H = p**2/2 + V(q) = p**2/2 - cos(q); V(x) = -cos(q)
    overlays a plot of V(q)''')
    h = 1e-3                # integrator's step size'
    if order == 3:
        map = Order3Map(h,dHdP,dHdX)
        text = u'Order3Map O(h\u00B3)'
    else:
        map = Order2Map(h,dHdP,dHdX)
        text = u'Order2Map O(h\u00B2)'
    q0   = radians(0)   # anfangswert q0
    # pmax=2
    pmax=1.9
    # proz = -1e-3          # increment pmax by proz%
    proz = 5.25          # increment pmax by proz%
    pmax=pmax*(1.+proz*1.e-2)
    pmin=pmax*1.e-3
    anz = 30           # anz trajectories
    dp = (pmax-pmin)/anz
    pvalues = [pmin+i*dp for i in range(anz,anz+1)]  # list of anfangswerte p0
    functions=[]
    for pv in pvalues:  # loop over several p
        # v = (pv,q0)
        v = (0.,q0,pv)
        function=[]
        for i in range(50000):   # loop over independent variable (time z.B.)
            pot = -cos(v[1])
            function.append((i,v[2]*70.,degrees(v[1]),pot))
            v = map.map(v)            # [p(t),q(t)] <---map---- [p(0),q(0)]
        functions.append(function)
    display(functions,text=text)
    print(
    '''
    =============================================================
    I understand now better the outcome of this excercise???!!!!
    =============================================================''')
def test2(fmap=11.3):              # use pendulum
    print(
    """"Shows Runge-Kutta integration of:
    Hamiltonian for pendulum
    H = p**2/2 + V(q) = p**2/2 - cos(q); V(q) = -cos(q)
    {p,-q} sind die kanonisch-konjugierten variablen
    q'' = - sin(q) ist die pendelgleichung; x=t die zeit
    """)
    h = 1e-3                # integrator's step size'
    if fmap == 11.2:
        map = RKMap11_2(h,pendel_f,pendel_g())
        text = 'RKMap11_2 O(h\u2075)'
    else:
        map = RKMap11_3(h,pendel_g())
        text = u'RKMap11_3 O(h\u2075)'
    q0   = radians(0)   # anfangswert q0
    pmax=1.9
    proz = 5.25      # increment pmax by proz%
    pmax=pmax*(1.+proz*1.e-2)
    pmin=pmax*1.e-3
    # anz = 10           # anz trajectories
    anz = 30           # anz trajectories
    dp = (pmax-pmin)/anz
    pvalues = [pmin+i*dp for i in range(anz,anz+1)]  # list of anfangswerte z0
    # print(pvalues)
    functions=[]
    for pv in pvalues:  # loop over several z
        v = (0.,q0,pv)
        function=[]
        for i in range(50000):   # loop over independent variable (time z.B.)
            pot = cos(v[2])
            function.append((i,v[2]*70.,degrees(v[1]),pot))
            v = map.map(v)            # [p(t),q(t)] <---map---- [p(0),q(0)]
        functions.append(function)
    display(functions,text=text)
    print(
    '''
    ===============================================================================
    I am still not 100% sure of correctness of the outcome of this excercise???!!!!
    ===============================================================================''')
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    test0()
    test1()
    test2()




