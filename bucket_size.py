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
from setup import CONF,Proton,objprnt,dictprnt
import matplotlib.pyplot as plt
from math import cos,pi,sqrt,sin,degrees,radians
from elements import RFG
import yaml
from fileLoader import unpack_list_of_dict

def display(functions):
    for function in functions:
        phi  = [x[0] for x in function]
        p1   = [x[1] for x in function]
        p2   = [x[2] for x in function]
        plt.plot(phi,p1,label='w(phi)',color='blue')
        plt.plot(phi,p2,label='',color='blue')
        plt.ylabel(r'$\Delta$w [Mev]')
        plt.xlabel(r'phi_sync [deg]')
    plt.title('longitudinal bucket')
    plt.show(block=True)

def psquared(H_invariant,phi,phis):
    ''' solves 0 = p**2 - V(phi) + H_invariant for p**2'''
    V = phi*cos(phis)-sin(phi)
    res = V-H_invariant
    return res

def read_yaml(filepath):          ## the principal YAML input parser
    with open(filepath,'r') as fileobject:
        in_data = yaml.load(fileobject)
    parameter_list = in_data['parameters']
    parameters     = unpack_list_of_dict(parameter_list)
#     print('parameters=\t',parameters)
    CONF['frequenz']         = parameters['frequency']
    CONF['injection_energy'] = parameters['TK_i']
    CONF['dP/P']             = parameters['dP/P'] * 1.e-2
    CONF['Ez_feld']          = parameters['Ez']
    CONF['soll_phase']       = parameters['phi_sync']
    CONF['dZ']               = parameters['dZ']
    CONF['spalt_laenge']     = parameters['gap']
    CONF['cavity_laenge']    = parameters['cav_len']
    CONF['wellenlänge']      = CONF['lichtgeschwindigkeit']/CONF['frequenz']
    CONF['spalt_spannung']   = CONF['Ez_feld']*CONF['spalt_laenge']
    CONF['input file']       = filepath
    return

def bucket(filepath):
    '''produce the longitudinal phase plots from Dr.Tiede'''
    read_yaml(filepath)

    phis=radians(CONF['soll_phase'])           # KNOB: soll phase

    dphi=1e-4                   # step size phase
    pmax=radians(+20.)          # phase upper limit
    pmin=radians(-40.)          # phase lower limit
    anz= int((pmax-pmin)/dphi)  # nboff phase steps

    if True:   # physics dimensions according to T.Wrangler pp.176
        ws=CONF['injection_energy']
        particle = Proton(ws)
        gapl=CONF['spalt_laenge']
        E0=CONF['Ez_feld']
        u0=CONF['spalt_spannung']
        fRF=CONF['frequenz']
        lamb=CONF['wellenlänge']
        rfg=RFG(U0=u0,PhiSoll=phis,fRF=fRF,label='RFG',gap=gapl,beam=particle,dWf=1.)
        dws=rfg.deltaW
        gammas=particle.gamma
        betas=particle.beta
        mc2=particle.e0
        q=1.
        T=particle.TrTf(gapl,fRF)
        A=2.*pi/(gammas*betas)**3/lamb
        B=q*E0*T/mc2
        p2w=sqrt(2.*B/A)*mc2   # conversion pk -> delta(w-ws) [Mev]

#         objprnt(rfg,text='cavity',filter=['matrix','beam'])
#         objprnt(particle,text='Particle')
        summary={
        '        cavity gap [m]':gapl,
        ' cavity frequency [Hz]':fRF,
        '  ref. energy Ws [MeV]':ws,
        '     delta-W/gap [MeV]':dws,
        '        particle gamma':gammas,
        '        particle  beta':betas,
        '         RF lambda [m]':lamb,
        '      cavity Ez [MV/m]':E0,
        '    m0*c**2 [MeV/c**2]':mc2,
        '                   ttf':T,
        '         particle type':particle.name,
        '            input file':CONF['input file']
        # '                     A':A,
        # '                     B':B,
        # '                   p2w':p2w,
        }
        dictprnt(summary,text='values for longitudinal dynamics')

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
    return

#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    import sys, os
    directory = os.path.dirname(__file__)
#     filepath = directory+'/fodo_with_10cav_per_RF.yml'       ## the default input file (YAML syntax)
    filepath = directory+'/fodo_with_10cav_per_RF(2).yml'       ## the default input file (YAML syntax)
    bucket(filepath)

