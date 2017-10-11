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
import sys
from math import sqrt
import numpy as np
from copy import copy

from setutil import PARAMS,DEBUG,Proton,tblprnt
from elements import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
import elements as ELM

## DEBUG MODULE
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF


class SIXD(ELM.D):
    """
    Sixtrack drift space
    """
    def __init__(self,length=0.,label="SIXD",viseo=0.,particle=PARAMS['sollteilchen'],position=[0.,0.,0.]):
        super().__init__(length=length, viseo=viseo, label='SIXD', particle=particle, position=position)
        self.off_soll = copy(self.particle)
    def shorten(self,l=0.):
        return SIXD(length=l,label=self.label, particle=self.particle, viseo=self.viseo)
    def adjust_energy(self,tkin):
        self.__init__(length=self.length, viseo=self.viseo, label=self.label, particle=self.particle(tkin), position=self.position)
        return self
    def map(self,i_track):
        return self._sixmap(i_track)
    def _sixmap(self,i_track):
        def fpsigma(psigma,soll):
            beta0 = soll.beta
            E0    = soll.e
            m0c2  = soll.e0
            res = (1+beta0**2*psigma)**2-(m0c2/E0)**2
            res = sqrt(res)/beta0-1.
            return res
        def einsplusfpsigma(psigma,soll):
            return 1.+fpsigma(psigma,soll)
        #conversion T3D ==> RipkenSchmidt (six)
        def t3d2six(i_track):
            soll     = self.particle
            x        = i_track[XKOO]       # [0]
            xp       = i_track[XPKOO]      # [1]
            y        = i_track[YKOO]       # [2]
            yp       = i_track[YPKOO]      # [3]
            z        = i_track[ZKOO]       # [4] z
            dp2p     = i_track[ZPKOO]      # [5] dp/p
            T        = i_track[EKOO]       # [6] summe aller delta-T
            s        = i_track[SKOO]       # [8] summe aller laengen

            E0       = soll.e
            beta0    = soll.beta
            p0       = soll.p          # cp-soll [MeV]
            m0c2     = soll.e0
            p        = p0/(1.-dp2p)
            E        = sqrt(p**2+m0c2**2) #E aus dp2p und p0
            tkin     = E-m0c2
            particle = self.off_soll(tkin=tkin)
            gb       = particle.gamma_beta
            beta     = particle.beta

            px       = gb*m0c2/E0*xp
            py       = gb*m0c2/E0*yp
            sigma    = z
            psigma   = ((beta0/beta/(1.-dp2p))-1.)/beta0**2
            f_track  = np.array([x,px,y,py,sigma,psigma,T,1.,s,1.])
            return f_track
        # conversion RipkenSchmidt (six) ==> T3D
        def six2t3d(i_track):
            soll     = self.particle
            x      = i_track[0]
            px     = i_track[1]
            y      = i_track[2]
            py     = i_track[3]
            sigma  = i_track[4]
            psigma = i_track[5]
            T      = i_track[EKOO]
            s      = i_track[SKOO]

            E0       = soll.e
            beta0    = soll.beta
            m0c2     = soll.e0
            eta      = beta0**2*psigma
            E        = (1.+eta)*E0
            tkin     = E-m0c2
            particle = self.off_soll(tkin=tkin)
            beta     = particle.beta
            gb       = particle.gamma_beta

            xp   = px/(gb*m0c2/E0)
            yp   = py/(gb*m0c2/E0)
            z    = sigma
            dp2p = 1.-beta0/beta/(1.+beta0**2*psigma)
            f_track = np.array([x,xp,y,yp,z,dp2p,T,1.,s,1.])
            return f_track
        # Ripken-Schnidt (six) map
        def rps_map(i_track,l):
            soll     = self.particle
            xi       = i_track[0]
            pxi      = i_track[1]
            yi       = i_track[2]
            pyi      = i_track[3]
            sigmai   = i_track[4]
            psigmai  = i_track[5]
            T        = i_track[EKOO]
            s        = i_track[SKOO]

            E0       = soll.e
            beta0    = soll.beta
            m0c2     = soll.e0
            eta      = beta0**2*psigmai
            E        = (1.+eta)*E0
            tkin     = E-m0c2
            particle = self.off_soll(tkin=tkin)
            beta     = particle.beta

            xf       = xi + pxi/einsplusfpsigma(psigmai,soll)*l
            pxf      = pxi
            yf       = yi + pyi/einsplusfpsigma(psigmai,soll)*l
            pyf      = pyi
            sigmaf  = sigmai + (1.-(beta0/beta)*(1.+0.5*(pxi**2+pyi**2)/einsplusfpsigma(psigmai,soll)**2))*l
            psigmaf = psigmai
            f_track = np.array([xf,pxf,yf,pyf,sigmaf,psigmaf,T,1.,s,1.])
            return f_track
        ##body
        f_track     = t3d2six(i_track)
        f_track     = rps_map(f_track,self.length)
        f_track     = six2t3d(f_track)
        f_track[SKOO] += self.length
        return f_track
    def soll_map(self,f_track):
        f_track[SKOO] += self.length
        return f_track

def test0():
    print('-----------------------------Test0---')        
    header_six=[
            'x',
            "px",
            'y',
            "py",
            'sigma',
            'psigma'
            ]
    header_t3d=[
            'x',
            "x'",
            'y',
            "y'",
            'z',
            'dp/p'
            ]

    tkin0      = 100.            #[MeV]
    l          = 50.e-3          #[m]
    sxdrift    = SIXD(length=l,particle=Proton(tkin=tkin0))
    
    # x,xp,y,yp,z,dp/p  T3D Kordinaten
    xi  = yi  = 1.e-3
    xpi = ypi = 1.e-3
    zi    = 1.e-3
    dp2pi = 1.e-2
    koord0 = np.array([xi,xpi,yi,ypi,zi,dp2pi,0.,1.,0.,1.])
    row = [[
            '{:9.6e}'.format(xi),
            '{:9.6e}'.format(xpi),
            '{:9.6e}'.format(yi),
            '{:9.6e}'.format(ypi),
            '{:9.6e}'.format(zi),
            '{:9.6e}'.format(dp2pi)
            ]]
    print('T3D IN\n'+tblprnt(header_t3d,row))
    
    koord = sxdrift.map(koord0)
    xf    = koord[0]
    xpf   = koord[1]
    yf    = koord[2]
    ypf   = koord[3]
    zf    = koord[4]
    dp2pf = koord[5]
    row = [[
            '{:9.6e}'.format(xf),
            '{:9.6e}'.format(xpf),
            '{:9.6e}'.format(yf),
            '{:9.6e}'.format(ypf),
            '{:9.6e}'.format(zf),
            '{:9.6e}'.format(dp2pf)
            ]]
    print('six OUT in T3D Koordinaten\n'+tblprnt(header_t3d,row))
##main
if __name__ == '__main__':
    test0()

