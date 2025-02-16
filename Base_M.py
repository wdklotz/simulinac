#!/Users/klotz/anaconda3/bin/python3.6
# -*- coding: utf-8 -*-
__version__='11.0.2.3'
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


import sys
import IGap
import unittest
import numpy   as NP
from copy import copy
from math import sqrt,degrees,cos
from setutil import I0,I1,WConverter,Proton
# from separatrix import w2phi
# from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, DSKOO, MDIM
# from setutil import DEBUG_ON,DEBUG_OFF,Proton
# from setutil import Ktp
# import warnings
# import math    as M
# import numpy   as NP
# import setutil as UTIL
# import OXAL as OX

twopi = 2.*pi

class Base_M(IGap.IGap):
    """ Base RF Gap-Model (A.Shishlo/J.Holmes ORNL/TM-2015/247) """
    def __init__(self):
        pass

    def configure(self,**kwargs):
        self.length       = 0. # 0. because it's a kick
        self.dWf          = FLAGS['dWf']
        self.mapping      = 'base'        # map model
        self.kwargs       = kwargs
        self.label        = 'BM' 

        self.particle  = kwargs.get('particle',Proton(50.))
        self.freq      = kwargs.get('freq',None)
        self.phisoll   = kwargs.get('phisoll',None)
        self.position  = kwargs.get('position',None)
        self.aperture  = kwargs.get('aperture',None)

        self.lamb      = PARAMS['clight']/self.freq
        self.ttf       = None
        self.qE0LT     = None
        self.deltaW    = None
        self.particlef = None  
        pass

    def values_at_exit(self):
        return dict(deltaw=self.deltaW,ttf=self.ttf,particlef=self.particlef,matrix=self.matrix)

    def map(self,i_track):
        return self.base_map(i_track)

    def toString(self):
        return 'base mapping: Base_M.base_map_1'

    def isAccelerating(self):
        return True

    def waccept(self,**kwargs): pass

    def register_mapper(self,master):
        master.register_mapping(self)
        pass

    def accept_register(self,master):
        self.master = master
        pass

    def adjust_energy(self, tkin):
        self.particle = Proton(tkin)
        self.OXAL_matrix()
        pass

    def base_map_1(self, i_track):
        """Neue map Version ab 03.02.2022 ist ein Remake um Korrecktheit der Rechnung zu testen. 
           Produziert dasselbe Verhalten wie base_map_0 """
        # def DEBUG_TRACK(inout,track):
        #     print('{} {} {}'.format('base_map',inout,track))
        # function body ================= function body ================= function body ================= 
        """ Mapping (i) to (O) in Base RF-Gap Model. (A.Shislo 4.2) """
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z
        zp       = i_track[ZPKOO]      # [5] dp/p
        T        = i_track[EKOO]       # [6] kinetic energy ref Teilchen
        S        = i_track[SKOO]       # [8] position gap

        particleRi = self.particle   # ref Teilchen (I)
        m0c2       = particleRi.e0
        betai      = particleRi.beta
        gammai     = particleRi.gamma
        gbi        = particleRi.gamma_beta
        wRi        = particleRi.tkin
        freq       = self.freq
        lamb       = self.lamb
        phisoll    = self.phisoll
        deg_phisoll= degrees(phisoll)
        qE0LT      = self.qE0LT
        deltaW     = self.deltaW
        
        # if 0: 
        #     DEBUG_ON()
        #     DEBUG_TRACK('tr_i',i_track)

        max_r  = 0.05              # max radial excursion [m]
        r      = sqrt(x**2+y**2)   # radial coordinate
        if r > max_r:
            raise OutOfRadialBoundEx(S)
        Kr     = (twopi*r)/(lamb*gbi)
        i0     = I0(Kr)            # bessel function I0
        i1     = I1(Kr)            # bessel function I1

        # if 0: print('Kr=',Kr,'r=',r,'gbi=',gbi,'i0=',i0,'i1=',i1)

        # ref Teilchen
        wRo = wRi + deltaW                           # ref Teilchen energy (O)
 
        # Teilchen
        converter   = WConverter(wRi,freq)
        deg_converter = degrees(converter.zToDphi(z)) 
        phiin       = converter.zToDphi(z) + phisoll 
        deg_phiin   = degrees(phiin)        # Teilchen phase (I)
        wo_wi       = qE0LT*i0*cos(phiin)                 # energy kick (Shislo 4.2.3)
        wi          = converter.Dp2pToDW(zp) + wRi        # Teilchen energy (I) dp/p --> dT
        wo          = wi + wo_wi                          # Teilchen energy (O)   
        dw          = wo - wRo                            # Differenz der energy kicks von Teilchen und ref Teilchen (entspricht delta**2)

        # DEBUG_OFF('base_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,i0,phis))

        particleRo = Proton(wRo)
        betao      = particleRo.beta
        gammao     = particleRo.gamma
        gbo        = particleRo.gamma_beta

        zo         = betao/betai*z                     # z (O) (4.2.5) A.Shishlo/J.Holmes
        zpo        = converter.DWToDp2p(dw)            # dW --> dp/p (O)

        # if 0: print('z ',z,'zpf ',zpf)

        factor = qE0LT/(m0c2*gbi*gbo)*i1               # common factor
        if r > 0.:
            xp  = gbi/gbo*xp - x/r*factor*M.sin(phiin)   # Formel 4.2.6 A.Shishlo/J.Holmes
            yp  = gbi/gbo*yp - y/r*factor*M.sin(phiin)
        elif r == 0.:
            xp  = gbi/gbo*xp
            yp  = gbi/gbo*yp

        f_track = NP.array([x, xp, y, yp, zo, zpo, T+deltaW, 1., S, 1.])

        # for DEBUGGING
        # if 0: DEBUG_TRACK('tr_f',f_track)
        # if 0:
        #     arrprnt([x*1.e3 for x in i_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
        #     arrprnt([x*1.e3 for x in f_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')

        # """ the parent reads these attributes below """
        self.particlef = particleRo
        return f_track
