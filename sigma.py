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
from math import pi,radians,degrees,sin,cos,sqrt
import numpy as NP
from copy import copy,deepcopy

from setutil import Proton, DEBUG, mxprnt, PARAMS, Ktw

DIM=6   # (0=x,1=x',2=y,3=y',4=z,5=dp/p) Trace3D

# DEBUG MODULE
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF

class Sigma(object):
    """ Utility class for handling the sigma-matrix """
    def __init__(self, twv0, epsx, epsy, epsz):
        self.matrix = NP.matrix(NP.zeros((DIM,DIM)))      ## sigma matrix (6x6)
        self._emitx = epsx
        self._emity = epsy
        self._emitz = epsz
        self._beam(twv0)
        
    def _beam(self,twv0):
        """ calc sigma-matrix from twiss-parameters """
        self.matrix[0,0] = self._emitx*twv0[Ktw.bx]
        self.matrix[2,2] = self._emity*twv0[Ktw.by]
        self.matrix[4,4] = self._emitz*twv0[Ktw.bz]

        self.matrix[1,1] = self._emitx*twv0[Ktw.gx]
        self.matrix[3,3] = self._emity*twv0[Ktw.gy]
        self.matrix[5,5] = self._emitz*twv0[Ktw.gz]
        
        self.matrix[0,1] = self.matrix[1,0] =  -self._emitx*twv0[Ktw.ax]
        self.matrix[2,3] = self.matrix[3,2] =  -self._emity*twv0[Ktw.ay]
        self.matrix[4,5] = self.matrix[5,4] =  -self._emitz*twv0[Ktw.az]
    
    def twiss(self):
        """ calc twiss-parameters from sigma-matrix """
        matrix = self.matrix
        bx  =  matrix[0,0]/self._emitx
        by  =  matrix[2,2]/self._emity
        bz  =  matrix[4,4]/self._emitz
        gx  =  matrix[1,1]/self._emitx
        gy  =  matrix[3,3]/self._emity
        gz  =  matrix[5,5]/self._emitz
        ax  = -matrix[0,1]/self._emitx
        ay  = -matrix[2,3]/self._emity
        az  = -matrix[4,5]/self._emitz
        return NP.array([bx,ax,gx,by,ay,gy,bz,az,gz])
    
    def sigv(self):
        """ calc envelopes from sigma-matrix """
        matrix = self.matrix
        sgx   = sqrt(matrix[0,0])   # sigmax  = <x*x>**1/2   [m]
        sgxp  = sqrt(matrix[1,1])   # sigmax' = <x'*x'>**1/2 [-]
        sgy   = sqrt(matrix[2,2])
        sgyp  = sqrt(matrix[3,3])
        sgz   = sqrt(matrix[4,4])
        sgzp  = sqrt(matrix[5,5])
        return NP.array([sgx,sgxp,sgy,sgyp,sgz,sgzp])

    def string(self):
        str = 'SIGMA:'
        for i in range(DIM):
            str += '\n'
            for k in range(DIM):
                str += '{:8.4g}  '.format(self.matrix[i,k])
        return str

    def RSRT(self,R):
        """
        Map this sigma-matrix through element R
        *) input R is ELM._matrix!
        *) returns the transformed Sigma object
        """
        # R(DIM,DIM) is same as R.__call__(DIM,DIM): upper (DIMxDIM) block matrix
        r6           = R(DIM,DIM)
        sigma        = deepcopy(self)       # IMPORTANT!!!
        sigma.matrix = r6 @ sigma.matrix @ r6.T # NP matrix multiplication
        return sigma

#todo: check eg_corr again - still with global delta-phi
    def apply_eg_corr(self,rf_gap, sigma_i, delta_phi, ksi=(0.,0.)):
        """
        Apply emmittance growth correction after passage through RF gap 
        ref: Appendix F Trace3D manual
        IN:
            self:      the sigma matrix (object) at gap exit a.k.a. sigma_f
            rf_gap:    the gap (object) for which the emiitance growth correction is added
            sigma_i:   the sigma matrix (object) at gap entrance
            delta_phi: the half-width of phase spread [rad]
            ksi:       tuple (x,y) as bunch center offset 
        OUT:
            self:      the corrected sigma matrix (object) at the gap exit a.k.a. sigma_f
        """
        def f(dp):
            res = ((sin(dp)-dp*cos(dp))/dp)
            res = res*3./(dp*dp)
            res = res - sin(dp)/dp
            res = (15.*res)/(dp*dp)
            return res
            
        def g(phis,dphi):
            res = 0.5*(1.+(sin(phis)**2-cos(phis)**2)*f(2.*dphi))
            return res
        
        Phis          = rf_gap.phis
        E0L           = rf_gap.EzAvg*rf_gap.gap
        T             = rf_gap.tr
        m0c2          = rf_gap.particle.e0
        lamb          = rf_gap.lamb
        particlei     = rf_gap.particle
        particlef     = rf_gap.particlef
        gamma_beta_f  = particlef.gamma_beta
        gamma_beta_av = (particlei.gamma_beta+gamma_beta_f)/2.
        kx            = -pi*E0L*T/(m0c2*gamma_beta_av**2*gamma_beta_f*lamb)
        cfactor1      = kx**2*(g(Phis,delta_phi)-(sin(Phis)*f(delta_phi))**2)
        ksix          = ksi[0]
        ksiy          = ksi[1]
        gamma_av      = (particlei.gamma+particlef.gamma)/2.
        kz            = -2.*kx*gamma_av**2*(1.+(delta_phi**2)/12.)
        cfactor2      = (kz*delta_phi)**2*((cos(Phis)**2)/8.+delta_phi*sin(Phis)/576.)
        delta_xp2_av  = cfactor1*(sigma_i.matrix[0,0]+ksix**2)
        delta_yp2_av  = cfactor1*(sigma_i.matrix[2,2]+ksiy**2)
        delta_dp2_av  = cfactor2*sigma_i.matrix[4,4]
        self.matrix[1,1] += delta_xp2_av
        self.matrix[3,3] += delta_yp2_av
        self.matrix[5,5] += delta_dp2_av
        DEBUG_MODULE('Phis ',degrees(Phis))
        DEBUG_MODULE('delta_phi ',degrees(delta_phi))
        DEBUG_MODULE('E0L ',E0L)
        DEBUG_MODULE('T ',T)
        DEBUG_MODULE('m0c2 ',m0c2)
        DEBUG_MODULE('lamb ',lamb)
        DEBUG_MODULE('gamma_beta_av ',gamma_beta_av)
        DEBUG_MODULE('gamma_beta_f ',gamma_beta_f)
        DEBUG_MODULE('gamma_av ',gamma_av)
        DEBUG_MODULE('kx ',kx)
        DEBUG_MODULE('kz ',kz)
        DEBUG_MODULE('cfactor1 ',cfactor1)
        DEBUG_MODULE('cfactor2 ',cfactor2)
        DEBUG_MODULE('ksi ',ksi)
        DEBUG_MODULE('delta_xp2_av ',delta_xp2_av)
        DEBUG_MODULE('delta_yp2_av ',delta_yp2_av)
        DEBUG_MODULE('delta_dp2_av ',delta_dp2_av)
        return self
def test0():
    print('-----------------------------Test0--')
    PARAMS['emitz_i'] = 0.0     # use this for test
    bx       = PARAMS['betax_i']
    ax       = PARAMS['alfax_i']
    gx       = (1+ax**2)/bx
    by       = PARAMS['betay_i']
    ay       = PARAMS['alfay_i']
    gy       = (1+ay**2)/bx
    twv0     = NP.array([bx,ax,gx,by,ay,gy])
    sg0      = Sigma(twv0)
    print(sg0.string())
def test1():
    from elements import RFG
    print('-----------------------------Test1--')

    particle = Proton(tkin=2.)
    R = RFG(particle=particle)

    bx       = PARAMS['betax_i']
    ax       = PARAMS['alfax_i']
    gx       = (1+ax**2)/bx
    by       = PARAMS['betay_i']
    ay       = PARAMS['alfay_i']
    gy       = (1+ay**2)/bx
    twv0     = NP.array([bx,ax,gx,by,ay,gy])
    sigma_i  = Sigma(twv0)
    s1 = sigma_i.matrix
    sigma_f = sigma_i.RSRt(R)   ## apply map to sigma
    s2 = sigma_f.matrix
    DEBUG('{sigma}\n',mxprnt(s1.A))
    DEBUG('{sigma_f} = {R}*{sigma}*{RT}\n',mxprnt(s2.A))

    sigma_fc = deepcopy(sigma_f)
    sigma_fc.apply_eg_corr(R,sigma_f,radians(25.),ksi=(0.01,0.01))
    s3 = sigma_fc.matrix
    DEBUG('{sigma_f}_corrected minus {sigma_f}_uncorrected\n',mxprnt((s3-s2).A))
# main ----------
if __name__ == '__main__':
    test0()    
    DEBUG_MODULE = DEBUG_ON
    test1()