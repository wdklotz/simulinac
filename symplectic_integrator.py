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
    def __init__(self,stepsize,phis):
        self.h = stepsize
        self.phis = phis
        self.cphis = cos(phis)
        self.c1 = 7./24.
        self.d1 = 2./3.
        self.c2 = 3./4.
        self.d2 = -2./3.
        self.c3 = -1./24.
        self.d3 = 1.
    def map(self,v0):
        v1 = self.step1(v0)
        v2 = self.step2(v1)
        v3 = self.step3(v2)
        return v3
    def g(self,p):
        return p
    def f(self,x):
        return cos(x)+self.cphis
    def step1(self,v):
        p0 = v[0]
        x0 = v[1]
        p1 = p0 + self.c1 * self.h * self.f(x0)
        x1 = x0 + self.d1 * self.h * self.g(p1)
        return (p1,x1)
    def step2(self,v):
        p1 = v[0]
        x1 = v[1]
        p2 = p1 + self.c2 *h * self.f(x1)
        x2 = x1 + self.d2 *h * self.g(p2)
        return (p2,x2)
    def step3(self,v):
        p2 = v[0]
        x2 = v[1]
        p3 = p2 + self.c3 * h * self.f(x2)
        x3 = x2 + self.d3 * h * self.g(p3)
        return (p3,x3)
def Vp(phi,phis):
    return +(sin(phi)-phi*cos(phis))
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
h = 1e-2
phis = radians(-90)   # parameter cos(phi0)
map = Order3Map(h,phis)
q0   = radians(90)   # anfangswert q(0)
pmax=+1.4134
anz = 30
dp = pmax/anz
pmin=-10*dp
pvalues = [pmin+i*dp for i in range(anz+1)]
functions=[]
for pv in pvalues:  # loop over several p
    v = (pv,q0)
    function=[]
    for i in range(int(10./h)):  # loop over independent varibale i.e. time
        pot = Vp(v[1],phis)
        function.append((i,v[0],degrees(v[1]),pot))
        v = map.map(v)  # [p(t),q(t)] <---map---- [p(0),q(0)]
    functions.append(function)
display(functions)
print(
'''
=============================================================
I still don't understand the outcome of this excercise???!!!!
=============================================================''')

    
    
    
