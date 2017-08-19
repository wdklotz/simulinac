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
from setutil import CONF,Proton
from matplotlib import pyplot as plt
from math import cos,pi,sqrt,sin,degrees,radians
from elements import RFG

def display_bucket(functions,phis,tki,gapl,qE0,fRF,name):
    # frame
    plt.figure('bucket size for '+CONF['lattice_version'],facecolor='#eaecef')
    # functions
    for function in functions:
        phi  = [x[0] for x in function]
        p1   = [x[1] for x in function]
        p2   = [x[2] for x in function]
        plt.plot(phi,p1,label='w(phi)',color='blue')
        plt.plot(phi,p2,label='',color='blue')
        plt.ylabel(r'$\Delta$W [KeV]')
        plt.xlabel(r'$\Delta$$\Phi$ [deg]')
    # make a text box
    xy_nx = plt.axis()
    tx1 = r'$\Phi$s: {:4.1f} [deg], '.format(phis)
    tx2 =    r'Tkin: {:4.1f} [MeV], '.format(tki)
    tx3 =       r'gap: {:4.1f} [mm] '.format(gapl*1.e3)
    tx4 =   r'Ez0: {:4.1f} [MeV/m], '.format(qE0)
    tx5 =     r'frf: {:4.1f} [MHz], '.format(fRF*1.e-6)
    tx6 =                      r'{}'.format(name)
    txt = tx1+tx2+tx3+'\n'+tx4+tx5+tx6
    plt.text(xy_nx[0]*0.75,xy_nx[3]*0.8,txt,bbox=dict(facecolor='bisque', alpha=0.8))
    # figure title
    plt.title('longitudinal bucket')
    plt.show(block=False)

def bucket():
    '''produce the longitudinal phase plots (Formeln T.Wangler pp.175)'''
    phis = radians(CONF['soll_phase'])           # KNOB: soll phase

    # Wertebereiche
    dphi = 1e-4                   # step size phase
    phimax = CONF['Dphimax+']     # stable phase upper limit
    phimin = CONF['Dphimax-']     # stable phase lower limit
    anz  =  int((phimax-phimin)/dphi) # nbof  phase steps

    tki      = CONF['injection_energy']
    particle = Proton(tkin=tki)
    gapl     = CONF['spalt_laenge']
    qE0      = CONF['Ez_feld']
    u0       = CONF['spalt_spannung']
    fRF      = CONF['frequenz']
    lamb     = CONF['wellenlänge']
    rfg      = RFG(U0=u0,PhiSoll=phis,fRF=fRF,label='RFG',gap=gapl,particle=particle,dWf=1.)
    gamma    = particle.gamma
    beta     = particle.beta
    m0c2     = particle.e0
    T        = particle.trtf(gapl,fRF)
    A        = 2.*pi/pow(gamma*beta,3)/lamb
    B        = qE0*T/m0c2
    H0       = -B*(sin(phis)-phis*cos(phis))      # hamiltonian
    H        = [H0 - (i-1)*H0/5. for i in range(11)]

    functions = []           # list of functions
    for h in H:
        function = []        # one function
        for i in range(-2,anz+2):
            phi = phimin+i*dphi
            w2 = 2./A*(h-B*(sin(phi)-phi*cos(phis)))
            if w2 < 0.:
                continue
            w = sqrt(w2)*m0c2*1.e3
            function.append([degrees(phi),w,-w])
        functions.append(function)
    display_bucket(functions,degrees(phis),tki,gapl,qE0,fRF,particle.name)
    return
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    print("bucket_size.py: sorry - nothing todo")
