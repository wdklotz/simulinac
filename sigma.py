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
from math import pi,radians,degrees,sin,cos
import numpy as np

from elements import RFG
from setutil import Proton,DEBUG,mxprnt

MDIM=6   #(0=x,1=x',2=y,3=y',4=z,5=dp/p) Trace3D

## DEBUG MODULE
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF

## sigma
class Sigma(object):
    def __init__(self,emitx=0.,betax=1.,alphax=0.,emity=0.,betay=1,alphay=0.,emitz=0.,betaz=1.,alphaz=0.):

        self.matrix = np.zeros((MDIM,MDIM))      ## sigma matrix (6x6)

        gammax = (1.+alphax**2)/betax
        gammay = (1.+alphay**2)/betay
        gammaz = (1.+alphaz**2)/betaz

        self.matrix[0,0] = emitx*betax
        self.matrix[2,2] = emity*betay
        self.matrix[4,4] = emitz*betaz

        self.matrix[1,1] = emitx*gammax
        self.matrix[3,3] = emity*gammay
        self.matrix[5,5] = emitz*gammaz
        
        self.matrix[0,1] = self.matrix[1,0] =  -emitx*alphax
        self.matrix[2,3] = self.matrix[3,2] =  -emity*alphay
        self.matrix[4,5] = self.matrix[5,4] =  -emitz*alphaz
    def string(self):
        str = 'SIGMA:'
        for i in range(MDIM):
            str += '\n'
            for k in range(MDIM):
                str += '{:8.4g}  '.format(self.matrix[i,k])
        return str
    def RSRt(self,R):
        """
        Map this sigma through element R
        """
        # remember R isinstance off ELM._matrix!
        new = Sigma()
        r = R(MDIM,MDIM) #!! use _matrix.__call__() to get element matrix as np.ndarray
        rt = np.transpose(r)
        s1 = self.matrix
        new.matrix = np.dot(r,np.dot(s1,rt))  ## matrix multiplication here
        return new
    def clone(self):
        """
        Clone new from sigma
        """
        new = Sigma()
        new.matrix = self.matrix.copy()
        return new
    def apply_eg_corr(self,rf_gap,sigma_i,delta_phi,ksi=(0.,0.)):
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
            # print('f()',res)
            return res
            
        def g(phis,dphi):
            res = 0.5*(1.+(sin(phis)**2-cos(phis)**2)*f(2.*dphi))
            # print('g()',res)
            return res
        
        Phis          = rf_gap.phis
        E0L           = rf_gap.u0
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
        # DEBUG_MODULE('Phis',degrees(Phis))
        # DEBUG_MODULE('delta_phi',degrees(delta_phi))
        # DEBUG_MODULE('E0L',E0L)
        # DEBUG_MODULE('T',T)
        # DEBUG_MODULE('m0c2',m0c2)
        # DEBUG_MODULE('lamb',lamb)
        # DEBUG_MODULE('gamma_beta_av',gamma_beta_av)
        # DEBUG_MODULE('gamma_beta_f',gamma_beta_f)
        # DEBUG_MODULE('gamma_av',gamma_av)
        # DEBUG_MODULE('kx',kx)
        # DEBUG_MODULE('kz',kz)
        # DEBUG_MODULE('cfactor1',cfactor1)
        # DEBUG_MODULE('cfactor2',cfactor2)
        # DEBUG_MODULE('ksi',ksi)
        # DEBUG_MODULE('delta_xp2_av',delta_xp2_av)
        # DEBUG_MODULE('delta_yp2_av',delta_yp2_av)
        # DEBUG_MODULE('delta_dp2_av',delta_dp2_av)
        return self
def test0():
    print('-----------------------------Test0--')
    print(Sigma(emitx=1.e-6,emity=1.e-3,emitz=1.).string())
    print(Sigma(emitx=1.e-6,alphax=1.,emity=1.e-3,alphay=1.,emitz=1.,alphaz=1.).string())
    print(Sigma(emitx=1.e-6,betax=10.,alphax=1.,emity=1.e-3,betay=10.,alphay=1.,emitz=1.,betaz=10.,alphaz=1.).string())
    s = Sigma()
    print('test __call__ method of SIGMA class:\n',s())
def test1():
    print('-----------------------------Test1--')

    particle = Proton(tkin=2.)
    R = RFG(particle=particle)

    sigma_i = Sigma(emitx=1.,betax=1.,alphax=1.,
                    emity=1.,betay=1.,alphay=0.,
                    emitz=1.,betaz=1.,alphaz=0.)
    s1 = sigma_i()
    sigma_f = sigma_i.RSRt(R)   ## apply map to sigma
    s2 = sigma_f()
    DEBUG('{sigma}\n',mxprnt(s1))
    DEBUG('{sigma_f} = {R}*{sigma}*{RT}\n',mxprnt(s2))

    sigma_fc = sigma_f.clone()
    sigma_fc.apply_eg_corr(R,sigma_f,radians(25.),ksi=(0.01,0.01))
    s3 = sigma_fc()
    DEBUG('{sigma_f}corrected minus {sigma_f}uncorrected\n',mxprnt(s3-s2))
## main ----------
if __name__ == '__main__':
    test0()    
    test1()