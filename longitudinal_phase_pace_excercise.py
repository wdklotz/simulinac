#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import Phys,Proton
from pylab import plot,show,legend
from math import cos,pi,sqrt,sin,fabs
'''produce the longitudinal phase plots from Dr.Tiede'''
def display(functions):
    for function in functions:
        phi  = [x[0] for x in function]
        p1   = [x[1] for x in function]
        p2   = [x[2] for x in function]
        plot(phi,p1,label='w(phi)',color='blue')
        plot(phi,p2,label='',color='blue')
    show()
def psquared(Hamiltonian,phi,phis):
    ''' solves 0 = p**2 - V(phi) + Hamiltonian for p**2'''
    V = phi*cos(phis)-sin(phi)
    res = V-Hamiltonian
    return res
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
rad=Phys['radians']
deg=Phys['degrees']
phis=-90.           # soll phase
phis=-60.

phis=rad*phis
dphi=1e-5           # step size phase
pmax=rad*91.        # phase upper limit
pmin=rad*-271.      # phase lower limit
anz= int((pmax-pmin)/dphi)  # nboff phase steps

functions=[]        # outer: list of functions
Hamilton=[-1.2+i*0.2 for i in range(10)]  
for HTW in Hamilton:
    function=[]        # inner: list of function values
    for i in range(anz):
        phi = pmin+i*dphi
        p = psquared(HTW,phi,phis)
        if p < 0.:
            continue
        p=sqrt(p)
        function.append([phi*deg,p,-p])
    functions.append(function)
display(functions)

    
    
    
