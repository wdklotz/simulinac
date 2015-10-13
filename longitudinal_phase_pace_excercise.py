#!/Users/klotz/pyzo2015a/python
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
from setup import CONF,Proton,objprnt,dictprnt
import matplotlib.pyplot as plt
from math import cos,pi,sqrt,sin,degrees,radians
from elements import RFG
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
phis=-30.
phis=radians(phis)

dphi=1e-4               # step size phase
# pmax=radians(91.)       # phase upper limit
# pmin=radians(-271.)     # phase lower limit
pmax=radians(30.)       # phase upper limit
pmin=radians(-80.)     # phase lower limit
anz= int((pmax-pmin)/dphi)  # nboff phase steps

with_physics_dimensions = True

if with_physics_dimensions:   # according to T.Wrangler pp.176
    # ws=CONF['injection_energy']
    ws=25.
    particle = Proton(ws)
    gapl=0.04
    u0=1./gapl
    fRF=1000.e6
    lamb=CONF['lichtgeschwindigkeit']/fRF
    rfg=RFG(U0=u0,PhiSoll=phis,fRF=fRF,label='RFG',gap=gapl,beam=particle,dWf=1.)
    dws=rfg.deltaW
    gammas=particle.gamma
    betas=particle.beta
    E0=u0*gapl
    mc2=particle.e0
    q=1.
    T=particle.TrTf(gapl,fRF)
    A=2.*pi/(gammas*betas)**3/lamb
    B=q*E0*T/mc2
    p2w=sqrt(2.*B/A)*mc2   # conversion pk -> delta(w-ws) [Mev]
    
    # objprnt(rfg,text='cavity',filter=['matrix','beam'])
    # objprnt(particle,text='Particle')
    summary={
    '        cavity gap [m]':gapl,
    ' cavity frequency [Hz]':fRF,
    '  ref. energy Ws [MeV]':ws,
    '     delta-W/gap [MeV]':dws,
    '        particle gamma':gammas,
    '        particle  beta':betas,
    '         RF lambda [m]':lamb,
    '      cavity Ez [MV/m]':E0,
    '   e-mc**2 {MeV/c**2]}':mc2,
    'time transition factor':T,
    '         particle type':particle.name,
    # '                     A':A,
    # '                     B':B,
    # '                   p2w':p2w,
    }
    dictprnt(summary,text='physics values')

# H_invariant=[-1.2+i*0.2 for i in range(7)]  
H_invariant=[-0.05+i*0.0025 for i in range(45)]  
    
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

    
    
    
