#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import Phys,Proton
import matplotlib.pyplot as plt
from math import cos,pi,sqrt,sin,degrees,radians
'''produce the longitudinal phase plots from Dr.Tiede'''
def display(functions):
    for function in functions:
        phi  = [x[0] for x in function]
        p1   = [x[1] for x in function]
        p2   = [x[2] for x in function]
        plt.plot(phi,p1,label='w(phi)',color='blue')
        plt.plot(phi,p2,label='',color='blue')
        plt.ylabel(r'$\Delta$w [Mev]')
        plt.xlabel(r'$\Delta$phi [deg]')
    plt.show()
def psquared(H_invariant,phi,phis):
    ''' solves 0 = p**2 - V(phi) + H_invariant for p**2'''
    V = phi*cos(phis)-sin(phi)
    res = V-H_invariant
    return res
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
phis=-90.           # KNOB: soll phase
phis=-60.
phis=radians(phis)

dphi=1e-4               # step size phase
pmax=radians(91.)       # phase upper limit
pmin=radians(-271.)     # phase lower limit
anz= int((pmax-pmin)/dphi)  # nboff phase steps

with_physics_dimensions = True

if with_physics_dimensions:   # according to T.Wrangler pp.176
    ws=Phys['injection_energy']
    # ws=25.
    particle = Proton(ws)
    gammas=particle.gamma
    betas=particle.beta
    lamb = Phys['wellenlänge']
    E0=Phys['spalt_spannung']
    mc2=particle.e0
    q=1.
    T=particle.TrTf()
    A=2.*pi/(gammas*betas)**3/lamb
    B=q*E0*T/mc2
    p2w=sqrt(2.*B/A)*mc2   # conversion pk -> delta(w-ws) [Mev]
    
    print('ws='.rjust(10),ws)
    print('gamma='.rjust(10),gammas)
    print('beta='.rjust(10),betas)
    print('lambda='.rjust(10),lamb)
    print('E0='.rjust(10),E0)
    print('mc2='.rjust(10),mc2)
    print('T='.rjust(10),T)
    print('A='.rjust(10),A)
    print('B='.rjust(10),B)
    print('p2w='.rjust(10),p2w)

H_invariant=[-1.2+i*0.2 for i in range(7)]  
    
functions=[]        # outer: list of functions
for HTW in H_invariant:
    function=[]        # inner: list of function values
    for i in range(anz):
        phi = pmin+i*dphi
        p = psquared(HTW,phi,phis)
        if p < 0.:
            continue
        p=p2w*sqrt(p)
        function.append([degrees(phi),p,-p])
    functions.append(function)
display(functions)

    
    
    
