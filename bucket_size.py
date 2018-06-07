#!/Users/klotz/anaconda3/bin/python3.6
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
from setutil import PARAMS,Proton
from matplotlib import pyplot as plt
from math import cos,pi,sqrt,sin,degrees,radians
from elements import RFG

def display_bucket(functions,phis,tkin,gap,EzAvg,freq,name):
    # frame
    # width=7; height=7
    # figsize = (width,height)
    # fighdr = 'bucket size for '+PARAMS['lattice_version']
    fig = plt.figure(num='bucket-size',facecolor='#eaecef')
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
    tx1 = r'$\Phi$s: {:4.1f} [deg], '. format(phis)
    tx2 =    r'Tkin: {:4.1f} [MeV], '. format(tkin)
    tx3 =       r'gap: {:4.1f} [mm]'.  format(gap*1.e3)
    tx4 =   r'Ez0: {:4.1f} [MeV/m], '. format(EzAvg)
    tx5 =     r'freq: {:4.1f} [MHz], '.format(freq*1.e-6)
    tx6 =                      r'{}'.  format(name)
    txt = tx1+tx2+tx3+'\n'+tx4+tx5+tx6
    plt.text(xy_nx[0]*0.75,xy_nx[3]*0.8,txt,bbox=dict(facecolor='bisque', alpha=0.8))
    # figure title
    plt.title('longitudinal bucket')

def bucket():
    '''produce the longitudinal phase plots (Formeln T.Wangler pp.175)'''
    phis = radians(PARAMS['phisoll'])           # KNOB: soll phase

    # Wertebereiche
    Dphi  = 1e-4                  # step size phase
    phi_2 = PARAMS['phi_2']       # stable phase lower limit
    psi   = PARAMS['psi']
    anz   = int(psi/Dphi)

    tkin     = PARAMS['injection_energy']
    particle = Proton(tkin=tkin)
    gap      = PARAMS['spalt_laenge']
    EzAvg    = PARAMS['EzAvg']
    freq     = PARAMS['frequenz']
    lamb     = PARAMS['wellenlänge']
    rfg      = RFG(EzAvg=EzAvg,PhiSoll=phis,fRF=freq,label='RFG',gap=gap,particle=particle,dWf=1.)
    gamma    = particle.gamma
    beta     = particle.beta
    m0c2     = particle.e0
    T        = particle.trtf(gap,freq)
    A        = 2.*pi/pow(gamma*beta,3)/lamb
    B        = EzAvg*T/m0c2
    H0       = -B*(sin(phis)-phis*cos(phis))      # hamiltonian
    H        = [H0 - (i-1)*H0/5. for i in range(11)]

    functions = []           # list of functions
    for h in H:
        function = []        # single function
        for i in range(-2,anz+2):
            phi = phi_2+i*Dphi
            w2 = 2./A*(h-B*(sin(phi)-phi*cos(phis)))
            if w2 < 0.:
                continue
            w = sqrt(w2)*m0c2*1.e3
            function.append([degrees(phi),w,-w])
        functions.append(function)
    display_bucket(functions,degrees(phis),tkin,gap,EzAvg,freq,particle.name)
    return
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    print("bucket_size.py: sorry - nothing todo")
