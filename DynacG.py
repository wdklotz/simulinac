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
#todo: revise slicing, mapping and soll-mapping
#todo: should not get negative kin. energy !!!!
import sys
import numpy as NP
from copy import copy
import math
from functools import partial
import warnings


from setutil import DEBUG, arrprnt, PARAMS, tblprnt
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO
from Ez0 import SFdata, Ipoly

# DEBUG__*
def DEBUG_ON(string,arg = '',end = '\n'):
    DEBUG(string,arg,end)
def DEBUG_OFF(string,arg = '',end = '\n'):
    pass
DEBUG_SLICE    = DEBUG_OFF
DEBUG_DYN_G    = DEBUG_OFF
DEBUG_TEST0    = DEBUG_ON

twopi = 2.*math.pi

""" Numerical computations in an accelerating gap or in a cavity;
    E.TANKE and S.VALERO
    3-Oct-2016 
"""
class _DYN_G(object):
    """ Wrapper to DYNAC's RF-gap model """
    def __init__(self, parent):
        def make_slices(parent, gap, SFdata, particle):
            """Slice the RF gap"""
            slices = []
            zl = -gap/2.*100.   # [m] --> [cm]
            zr = -zl
            E0z = 0.
            z = 0.
            for poly in SFdata.Ez_poly:
                zil = poly.zl
                zir = poly.zr
                if zil < zl or zir > zr: continue
                # instanciate _DYN_Gslices
                slice = _DYN_Gslice(parent, poly, particle)
                slices.append(slice)
            return slices

        def configure_slices(slices, phis, tkin):
            next_phase = phis
            next_tkin  = tkin
            Tklist = []
            for slice in slices:
                print(next_tkin)
                # setting phase and energy @ slice entrance
                slice.adjust_slice_parameters(next_phase, next_tkin)
                Tklist.append(slice.Tk)
                next_phase = slice.PHOUT
                next_tkin  = slice.WOUT
            deltaW  = next_tkin-tkin        # total energy kick as sum over slices
            tr      = NP.sum(NP.array(Tklist))/len(Tklist)
            return deltaW,tr

        # _DYN_G attributes
        self.phis     = parent.phis
        self.freq     = parent.freq
        self.gap      = parent.gap
        self.dWf      = parent.dWf
        self.lamb     = parent.lamb
        self.SFdata   = parent.SFdata
        self.matrix   = parent.matrix
        self.particle = parent.particle
        if self.SFdata == None:
            raise RuntimeError('_DYN_G: missing E(z) table - STOP!')
            sys.exit(1)
        else:
            # slice the gap
            self.slices = \
                make_slices(self, self.gap, self.SFdata, self.particle)
            DEBUG_DYN_G('_DYN_G:make_slices()\n',list(self.slices))
            # configure slices
            self.deltaW, self.tr = \
                configure_slices(self.slices, self.phis, self.particle.tkin)
            # UPDATE linear matrix with this deltaW
            self.matrix[EKOO, DEKOO] = self.deltaW

    def map(self,i_track):
        """ Mapping from (i) to (f) """
        for slice in self.slices:
            f_track = slice.slice_map(i_track)   # map slice with DYNAC gap-model
            i_track = f_track
        
        # parent property
        self.particlef = copy(self.particle)(self.particle.tkin + self.deltaW)

        # for DEBUGGING
        if DEBUG_DYN_G == DEBUG_ON:
            itr = i_track.copy()
            ftr = f_track.copy()
            for i in range(len(f_track)-4):
                itr[i]  = itr[i]*1.e3
                ftr[i]  = ftr[i]*1.e3
            arrprnt(itr, fmt = '{:6.4g},', txt = 'dynac_map:i_track:')
            arrprnt(ftr, fmt = '{:6.4g},', txt = 'dynac_map:f_track:')

        return f_track
        
    def soll_map(self, i_track):
        return self.map(i_track)

class _DYN_Gslice(object):
    """ The DYNAC gap-model """
    def __init__(self, parent, poly, particle):
        # _DYN_Gslice attributes (inherited)
        self.parent   = parent
        self.freq     = parent.freq     # frquency
        self.omega    = twopi*self.freq # Kreisfrequenz
        self.lamb     = parent.lamb     # Wellenlaenge
        self.particle = copy(particle)  # soll particle (copy)
        self.poly     = poly            # the current interval
        self.SFdata   = parent.SFdata   # reference to superfish data
    #todo: calculate Tk better ??
        self.Tk       = 0.85

    def time_array(self,betac,h,zarr):
        """ create arrival times azimutal positions 
            IN-->betac = velocity
                 h = azimutal step size
                 zarr = azimutal positions
            
            OUT<--array of arrival times
        """
        th = - h / betac         # [s] time step size
        # t0 is time at slice entry!
        t0 = - zarr[0] / betac   # [s] time on interval entry
        t4 = th+t0
        t3 = 0.75*th+t0
        t2 = 0.50*th+t0
        t1 = 0.25*th+t0
        return NP.array([t0,t1,t2,t3,t4])
        
    def adjust_slice_parameters(self, phin, tkin):
        """ Adjust soll-energy-dependent parameters for this slice 
            IN-->phin = phase of soll at entry
                 tkin = kin.energy at entry
        """
        # new soll-particle energy
        particle   = self.particle(tkin)
        m0c3       = particle.m0c3
        betac      = particle.betac

        self.h = 1.e-2*2.*self.poly.dz # [m] azimutal step size
        # z0 = z at slice entry!
        h  = self.h
        z0 = 1.e-2*self.poly.zl        # [m] interval left border
        z4 = h+z0                      # [m] interval right border
        z3 = 0.75*h+z0
        z2 = 0.50*h+z0
        z1 = 0.25*h+z0
        self.zarr = NP.array([z0,z1,z2,z3,z4])                # [m]
        tarr      = self.time_array(betac, h, self.zarr)      # [s]
        I1        = self.Integral1(self.zarr, tarr, h, self.omega, phin, particle)
        I3        = self.Integral3(self.zarr, tarr, h, self.omega, phin, particle)
        # soll parameter
        self.BETA     = particle.beta                     # beta
        self.GAMMA    = particle.gamma                    # gamma
        self.WIN      = tkin                              # (z0) kin.energy
        self.PHIN     = phin                              # (z0) phase
        self.deltaW   = I1                                # (z4) delta-kin.energy
        self.deltaT   = I3/m0c3                           # (z4) delta-time
        self.deltaP   = self.deltaT*self.omega            # (z4) delta-phase
        self.WOUT     = tkin + self.deltaW                # (z4) kin.energy
        self.PHOUT    = phin + self.deltaP                # (z4) phase

        DEBUG_SLICE('_DYN_Gslice: {}\n'.format(self),self.__dict__)
        DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():zarr.........[m]: ', self.zarr)
        DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():tarr......[psec]: ', 1.e12*tarr)
        # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():K1......[1/m**2]: ', K1)
        DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():I1..........[MV]: ', I1)
        # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():I2........[MV*m]: ', self.I2)
        DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():I3........[MV*m]: ', I3)
        # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():I4.....[MV*m**2]: ', self.I4)
        # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():J1.........[1/m]: ', self.J1)
        # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():J2............[]: ', self.J2)
        # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():J3...........[m]: ', self.J3)
        DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():tkin.......[MeV]: ', tkin)
        DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():phin..........[]: ', phin)
        DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():WOUT.......[MeV]: ', self.WOUT)
        DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():PHOUT.........[]: ', self.PHOUT)
        DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():deltaP........[]: ', self.deltaP)
        # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():deltaT.....[sec]: ', self.deltaT)
        DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():deltaW.....[Kev]: ', self.deltaW*1.e3)
        DEBUG_SLICE('============================================================ adjust_slice_parameters() end')
        return

    def slice_map(self, i_track):
        """ The DYNAC RF-gap model map of a slice """
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z
        zp       = i_track[ZPKOO]      # [5] dp/p
        T        = i_track[EKOO]       # [6] summe aller dT
        S        = i_track[SKOO]       # [8] summe aller laengen

        # particle kin.energy and phase at entrance (z0)
        win      = (zp * (self.GAMMA+1.)/self.GAMMA +1.) * self.WIN      # (z0) W
        pin      = - z*twopi/(self.BETA*self.lamb) + self.PHIN           # (z0) phase
        particle = copy(self.particle)(win)                              # (z0) particle

        # aliases
        gamma   = particle.gamma
        beta    = particle.beta
        gambeta = particle.gamma_beta
        lamb    = self.lamb
        m0c2    = particle.m0c2
        m0c3    = particle.m0c3
        betac   = particle.betac
        omega   = self.omega
        gbroot  = math.sqrt(gambeta)
        gb3     = gambeta**3
        K1      = (math.pi/self.lamb)**2/gb3  # [1./m**2] common factor
 
        # Integrale
        tarr = self.time_array(betac, self.h, self.zarr)
        I1   = self.Integral1(self.zarr, tarr, self.h, omega, pin, particle)
        I2   = self.Integral2(self.zarr, tarr, self.h, omega, pin, particle)
        I3   = self.Integral3(self.zarr, tarr, self.h, omega, pin, particle)
        I4   = self.Integral4(self.zarr, tarr, self.h, omega, pin, particle)
        J1   = self.Jntegral1(self.zarr, tarr, self.h, omega, pin, particle)
        J2   = self.Jntegral2(self.zarr, tarr, self.h, omega, pin, particle)
        J3   = self.Jntegral3(self.zarr, tarr, self.h, omega, pin, particle)

        # Picht transformation
        r   = NP.array([x,y])                      # =(x,y)
        rp  = NP.array([xp,yp])                    # =(x',y')
        R   = r*gbroot                             # =(X,Y)
        Rp  = rp*gbroot+0.5*R*gamma/(gamma**2-1.)  # =(X',Y')
        # delta gamma
        dgamma  = ((1. + NP.dot(R,R)*K1)*I1 + NP.dot(R,Rp)*K1*I2)/m0c2
        deltaW  = dgamma * m0c2                    # (z4) delta-kin.energy
        # delta time
        dtime  = ((1. + NP.dot(R,R)*K1)*I3 + NP.dot(R,Rp)*K1*I4)/m0c3
        dphi   = dtime * omega                     # (z4) delta-phase
        dz1    = - beta*lamb/twopi*dphi            # z = distance from soll
        dz2    = - betac * dtime                   # z = distance from soll
        if abs(dz1-dz2) > 1.e-15:
            warnings.warn('|delta-z difference| too large - should be equal')
        DEBUG_OFF('(deltaW[KeV], dphi[mdeg]) ',(deltaW*1.e3, math.degrees(dphi)*1.e3))

        # mapping of reduced coordinates
        dR     = R*J2 + Rp*J3  # delta-radius
        dRp    = R*J1 + Rp*J2  # delta-radius primed  
        Rf     = R + dR
        Rpf    = Rp + dRp

        # Picht back-transformation
        particle = copy(self.particle)(tkin = win + deltaW)
        gamma    = particle.gamma
        gambeta  = particle.gamma_beta
        gbroot   = math.sqrt(gambeta)
        rf       = Rf/gbroot
        rpf      = (Rpf - 0.5*Rf*gamma/(gamma**2-1.)) / gbroot

        # new z
        pout     = pin + dphi
        dp       = pout - self.PHOUT
        zf       = -beta*self.lamb/twopi*dp     # z out
        # new dp/p 
        wout     = win + deltaW
        dw       = wout - self.WOUT
        zpf      = gamma/(gamma+1.)*dw/wout     # delta-p/p out
  
        # new track point
        f_track = NP.array([x,xp,y,yp,z,zp,T,1.,S,1.])
        f_track[XKOO]  = rf[0]
        f_track[YKOO]  = rf[1]
        f_track[XPKOO] = rpf[0]
        f_track[YPKOO] = rpf[1]
        f_track[ZKOO]  = zf
        f_track[ZPKOO] = zpf
        f_track[EKOO]  = T + deltaW

        # DEBUG_SLICE('_DYN_Gslice:slice_map():R............[m]: ', R)
        # DEBUG_SLICE('_DYN_Gslice:slice_map():Rp.........[rad]: ', Rp)
        # DEBUG_SLICE('_DYN_Gslice:slice_map():dR..........[um]: ', dR*1e6)
        # DEBUG_SLICE('_DYN_Gslice:slice_map():dRp.......[urad]: ', dRp*1e6)
        # DEBUG_SLICE('_DYN_Gslice:slice_map():dgamma..........: ', dgamma)
        # DEBUG_SLICE('_DYN_Gslice:slice_map():dtime.....[psec]: ', dtime*1.e12)
        # DEBUG_SLICE('_DYN_Gslice:slice_map():r=(x,y)......[m]: ', rf)
        # DEBUG_SLICE("_DYN_Gslice:slice_map():rp=(x',y').[rad]: ", rpf)
        # DEBUG_SLICE('============================================================ slice_map() end')
        return f_track

    def Integral1(self, zarr, tarr, h, omega, phin, particle):
        coeff = NP.array([7., 32., 12., 32., 7.])
        res   = 0.
        # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
        E = partial(self.SFdata.Ez0t, omega = omega, phi = phin)
        for i in range(0,len(coeff)):
            z = 1.e2*zarr[i]        # [m] --> [cm]
            res = res + coeff[i] * E(z, tarr[i])
            # DEBUGING
            # if i == 0: print('====E(z,0,t)')
            # t = tarr[i]
            # cos = math.cos(omega*t+phin)
            # Ez = Ipoly(z,self.SFdata.Ez_poly)
            # probe = E(z,tarr[i]) - Ez*cos
            # print('(z,t,I1,E,cos,0.0) ',(z,t,res,Ez,cos,probe))
        res = res * h / 90.
        return res

    def Integral2(self, zarr, tarr, h, omega, phin, particle):
        coeff = NP.array([0., 8., 6., 24., 7.])
        res   = 0.
        # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
        E = partial(self.SFdata.Ez0t, omega = omega, phi = phin)
        for i in range(1,len(coeff)):
            z = 1.e2*zarr[i]        # [m] --> [cm]
            res = res + coeff[i] * E(z, tarr[i])
        res = res * h**2 / 90.
        return res

    def Integral3(self, zarr, tarr, h, omega, phin, particle):
        coeff = NP.array([0., 8., 6., 24., 7.])
        res   = 0.
        bg3 = particle.gamma_beta**3
        # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
        E = partial(self.SFdata.Ez0t, omega = omega, phi = phin)
        for i in range(1,len(coeff)):
            z = 1.e2*zarr[i]        # [m] --> [cm]
    #todo: use beta_gamma[i], i=1,4,1
            res = res + coeff[i]/bg3 * E(z, tarr[i])
        res = res * h**2 / 90.
        return res

    def Integral4(self, zarr, tarr, h, omega, phin, particle):
        coeff = NP.array([0., 2., 3., 18., 7.])
        res   = 0.
        bg3 = particle.gamma_beta**3
        # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
        E = partial(self.SFdata.Ez0t, omega = omega, phi = phin)
        for i in range(1,len(coeff)):
            z = 1.e2*zarr[i]        # [m] --> [cm]
    #todo: use beta_gamma[i], i=1,4,1
            res = res + coeff[i]/bg3 * E(z, tarr[i])
        res = res * h**3 / 90.
        return res

    def Jntegral1(self, zarr, tarr, h, omega, phin, particle):
        coeff = NP.array([7., 32., 12., 32., 7.])
        res   = 0.
        g2m1 = particle.gamma**2-1.
    #todo: use G1(gamma[i]), i=1,4,1
        G1    = g2m1**(-1.5)/(2.*particle.m0c3)
        # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
        Ep = partial(self.SFdata.dEz0tdt, omega = omega, phi = phin)
        # print('====dE(z,0,t)/dt')
        for i in range(0,len(coeff)):
            z = 1.e2*zarr[i]        # [m] --> [cm]
            res = res + coeff[i] * G1 * Ep(z, tarr[i])
            # print('(z,Ep(z, tarr[i])) ',(z,Ep(z, tarr[i])))
        res = res * h / 90.
        return res

    def Jntegral2(self, zarr, tarr, h, omega, phin, particle):
        coeff = NP.array([0., 8., 6., 24., 7.])
        res   = 0.
        g2m1 = particle.gamma**2-1.
    #todo: use G1(gamma[i]), i=1,4,1
        G1    = g2m1**(-1.5)/(2.*particle.m0c3)
        # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
        Ep = partial(self.SFdata.dEz0tdt, omega = omega, phi = phin)
        for i in range(1,len(coeff)):
            z = 1.e2*zarr[i]        # [m] --> [cm]
            res = res + coeff[i] * G1 * Ep(z, tarr[i])
        res = res * h**2 / 90.
        return res

    def Jntegral3(self, zarr, tarr, h, omega, phin, particle):
        coeff = NP.array([0., 2., 3., 18., 7.])
        res   = 0.
        g2m1 = particle.gamma**2-1.
    #todo: use G1(gamma[i]), i=1,4,1
        G1    = g2m1**(-1.5)/(2.*particle.m0c3)
        # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
        Ep = partial(self.SFdata.dEz0tdt, omega = omega, phi = phin)
        for i in range(1,len(coeff)):
            z = 1.e2*zarr[i]        # [m] --> [cm]
            res = res + coeff[i] * G1 * Ep(z, tarr[i])
        res = res * h**3 / 90.
        return res

def test0():
    import elements as ELM
    print('-----------------------------------TEST 0----------------')
    print('test _DYN_Gslice:slice_map()...')
    input_file='SF_WDK2g44.TBL'
    Ezpeak = 1.4
    SF_tab = SFdata(input_file,Ezpeak)
    x  = 1.0e-2
    xp = 1.0e-3
    y  = 1.0e-2
    yp = 1.0e-3
    z  = 0.0
    zp = 0.0
    i_track = NP.array([ x, xp, y, yp, z, zp, PARAMS['sollteilchen'].tkin, 1., 0., 1.])

    dyng = ELM.RFG(gap = 0.048,SFdata = SF_tab,mapping = 'dyn')
    DEBUG_TEST0('_DYN_Gslice:test0():i_track:\n', str(i_track))
    f_track = dyng.map(i_track)
    DEBUG_TEST0('_DYN_Gslice:test0():f_track:\n', str(f_track))

if __name__ == '__main__':
    test0()
