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
import unittest

from setutil import Proton, mxprnt, PARAMS, Ktw, DEBUG_ON, DEBUG_OFF

DIM = 6   # (0=x,1=x',2=y,3=y',4=z,5=dp/p) Trace-3D

class Sigma(object):
    """ class for the sigma-matrix """
    def __init__(self, twiss_vec0, epsx, epsy, epsz):
        self.matrix = NP.zeros((DIM,DIM))      ## sigma matrix (6x6)
        self.emitx = epsx
        self.emity = epsy
        self.emitz = epsz
        """ sigma-matrix from initial twiss-parameters """
        self.matrix[0,0] = self.emitx*twiss_vec0[Ktw.bx]
        self.matrix[2,2] = self.emity*twiss_vec0[Ktw.by]
        self.matrix[4,4] = self.emitz*twiss_vec0[Ktw.bz]

        self.matrix[1,1] = self.emitx*twiss_vec0[Ktw.gx]
        self.matrix[3,3] = self.emity*twiss_vec0[Ktw.gy]
        self.matrix[5,5] = self.emitz*twiss_vec0[Ktw.gz]
        
        self.matrix[0,1] = self.matrix[1,0] =  -self.emitx*twiss_vec0[Ktw.ax]
        self.matrix[2,3] = self.matrix[3,2] =  -self.emity*twiss_vec0[Ktw.ay]
        self.matrix[4,5] = self.matrix[5,4] =  -self.emitz*twiss_vec0[Ktw.az]
    def sig_twiss_vec_get(self):
        """ return twiss-vector from sigma-matrix """
        bx  =  self.matrix[0,0]/self.emitx
        by  =  self.matrix[2,2]/self.emity
        bz  =  self.matrix[4,4]/self.emitz
        gx  =  self.matrix[1,1]/self.emitx
        gy  =  self.matrix[3,3]/self.emity
        gz  =  self.matrix[5,5]/self.emitz
        ax  = -self.matrix[0,1]/self.emitx
        ay  = -self.matrix[2,3]/self.emity
        az  = -self.matrix[4,5]/self.emitz
        return NP.array([bx,ax,gx,by,ay,gy,bz,az,gz])
    def sig_sigma_vec_get(self):
        """ return sigma-vector from sigma-matrix """
        sgx   = sqrt(self.matrix[0,0])   # sigmax  = <x*x>**1/2   [m]
        sgxp  = sqrt(self.matrix[1,1])   # sigmax' = <x'*x'>**1/2 [-]
        sgy   = sqrt(self.matrix[2,2])
        sgyp  = sqrt(self.matrix[3,3])
        sgz   = sqrt(self.matrix[4,4])
        sgzp  = sqrt(self.matrix[5,5])
        return NP.array([sgx,sgxp,sgy,sgyp,sgz,sgzp])

def sig_apply_eg_corr(rf_gap, sigma_i, delta_phi, ksi=(0.,0.)):
    """
    Apply emmittance growth correction after passage through RF gap 
    ref: Appendix F Trace3D manual
    IN:
        rf_gap:    the gap object for which the correction is applied
        sigma_i:   the sigma object at gap entrance
        delta_phi: the half-width of phase spread [rad]
        ksi:       tuple (x,y) as bunch center offset 
    OUT:
        sigma_f:   the corrected sigma object at the gap exit
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
    
    Phis          = rf_gap.phisoll
    E0L           = rf_gap.EzAvg*rf_gap.gap
    ttf           = rf_gap.ttf
    m0c2          = rf_gap.particle.e0
    lamb          = rf_gap.lamb
    particlei     = rf_gap.particle
    particlef     = rf_gap.particlef
    gamma_beta_f  = particlef.gamma_beta
    gamma_beta_av = (particlei.gamma_beta+gamma_beta_f)/2.
    kx            = -pi*E0L*ttf/(m0c2*gamma_beta_av**2*gamma_beta_f*lamb)
    cfactor1      = kx**2*(g(Phis,delta_phi)-(sin(Phis)*f(delta_phi))**2)
    ksix          = ksi[0]
    ksiy          = ksi[1]
    gamma_av      = (particlei.gamma+particlef.gamma)/2.
    kz            = -2.*kx*gamma_av**2*(1.+(delta_phi**2)/12.)
    cfactor2      = (kz*delta_phi)**2*((cos(Phis)**2)/8.+delta_phi*sin(Phis)/576.)
    delta_xp2_av  = cfactor1*(sigma_i.matrix[0,0]+ksix**2)
    delta_yp2_av  = cfactor1*(sigma_i.matrix[2,2]+ksiy**2)
    delta_dp2_av  = cfactor2*sigma_i.matrix[4,4]
    
    sigma         = deepcopy(sigma_i)     # the new sigma object
    sigma.matrix[1,1] += delta_xp2_av
    sigma.matrix[3,3] += delta_yp2_av
    sigma.matrix[5,5] += delta_dp2_av
    DEBUG_OFF('Phis {}'.format(degrees(Phis)))
    DEBUG_OFF('delta_phi {}'.format(degrees(delta_phi)))
    DEBUG_OFF('E0L {}'.format(E0L))
    DEBUG_OFF('ttf {}'.format(ttf))
    DEBUG_OFF('m0c2 {}'.format(m0c2))
    DEBUG_OFF('lamb {}'.format(lamb))
    DEBUG_OFF('gamma_beta_av {}'.format(gamma_beta_av))
    DEBUG_OFF('gamma_beta_f {}'.format(gamma_beta_f))
    DEBUG_OFF('gamma_av {}'.format(gamma_av))
    DEBUG_OFF('kx {}'.format(kx))
    DEBUG_OFF('kz {}'.format(kz))
    DEBUG_OFF('cfactor1 {}'.format(cfactor1))
    DEBUG_OFF('cfactor2 {}'.format(cfactor2))
    DEBUG_OFF('ksi {}'.format(ksi))
    DEBUG_OFF('delta_xp2_av {}'.format(delta_xp2_av))
    DEBUG_OFF('delta_yp2_av {}'.format(delta_yp2_av))
    DEBUG_OFF('delta_dp2_av {}'.format(delta_dp2_av))
    return sigma
def sig_map(sigma_i,R):
    """
    Map this sigma-matrix through element R
    *) input R is ELM._matrix!
    *) returns the transformed Sigma object
    """
    # R(DIM,DIM) is same as R.__call__(DIM,DIM): upper (DIMxDIM) block matrix
    r6             = R(DIM,DIM)
    sigma_f        = deepcopy(sigma_i)          # the new sigma object
    sigma_f.matrix = r6 @ sigma_i.matrix @ r6.T # NP matrix multiplication
    return sigma_f
    
class TestElementMethods(unittest.TestCase):
    def test_instanciate_sigma(self):
        print('---------------------test_instanciate_sigma--')
        PARAMS['emitz_i'] = 0.0     # use this for test
        bx       = 1.
        ax       = 0.
        gx       = (1+ax**2)/bx
        by       = 1.
        ay       = 0.
        gy       = (1+ay**2)/bx
        bz       = 1.
        az       = 0.
        gz       = (1+ay**2)/bx
        twiss_vec0 = NP.array([bx,ax,gx,by,ay,gy,bz,az,gz])
        sg0      = Sigma(twiss_vec0,1.,1.,1.)
        if DEBUG_ON('[sigma]'):
            print('{}'.format(mxprnt(sg0.matrix)))
    def test_map_sigma(self):
        print('---------------------test_map_sigma--')
        from elements import RFG
        particle = Proton(2.)
        R = RFG('test-gap',5.0, -30.,0.022, 816.E6, particle=particle)
        bx       = 1.
        ax       = 0.
        gx       = (1+ax**2)/bx
        by       = 1.
        ay       = 0.
        gy       = (1+ay**2)/bx
        bz       = 1.
        az       = 0.
        gz       = (1+ay**2)/bx
        twiss_vec0 = NP.array([bx,ax,gx,by,ay,gy,bz,az,gz])
        sigma_i    = Sigma(twiss_vec0,1.,1.,1.)
        simx       = sigma_i.matrix
        sigma_f    = sig_map(sigma_i,R)   ## apply map to sigma
        sfmx       = sigma_f.matrix
        if DEBUG_ON('[sigma_i]'):
            print('{}'.format(mxprnt(simx)))
        if DEBUG_ON('[sigma_f] = [R]*[sigma]*[R]^T'):
            print('{}'.format(mxprnt(sfmx)))
    def test_eg_correction(self):
        print('---------------------test_eg_correction--')
        from elements import RFG
        particle = Proton(2.)
        R = RFG('test-gap',5.0, -30.,0.022, 816.E6, particle=particle)
        bx       = 1.
        ax       = 0.
        gx       = (1+ax**2)/bx
        by       = 1.
        ay       = 0.
        gy       = (1+ay**2)/bx
        bz       = 1.
        az       = 0.
        gz       = (1+ay**2)/bx
        twiss_vec0 = NP.array([bx,ax,gx,by,ay,gy,bz,az,gz])
        sigma_i    = Sigma(twiss_vec0,1.,1.,1.)
        sigma_f    = sig_map(sigma_i,R)   ## map to sigma through R
        sfmx       = sigma_f.matrix
        delta_phi  = radians(5.)
        sigma_fc = sig_apply_eg_corr(R,sigma_f,delta_phi,ksi=(0.01,0.01))
        sfcmx = sigma_fc.matrix
        if DEBUG_ON('[sigma_f]-corrected minus [sigma_f]-uncorrected'):
            print('{}'.format(mxprnt((sfcmx-sfmx))))
# main ----------
if __name__ == '__main__':
    unittest.main()
