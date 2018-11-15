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
import math as MATH
from functools import partial
import warnings

from setutil import DEBUG, arrprnt, PARAMS, tblprnt, Ktp, WConverter
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO
from Ez0 import SFdata, Ipoly, DynacSteps

# DEBUG__*
def DEBUG_ON(string,arg = '',end = '\n'):
    DEBUG(string,arg,end)
def DEBUG_OFF(string,arg = '',end = '\n'):
    pass
# DEBUG_SLICE    = DEBUG_OFF
DEBUG_DYN_G    = DEBUG_OFF

twopi = 2.*MATH.pi

""" Numerical computations in an accelerating gap or in a cavity;
    E.TANKE and S.VALERO
    3-Oct-2016 
"""

class StepFactory(object):
    """ 
    StepFactory
    """
    def __init__(self,gap,SFdata):
        self.zl = -gap*100./2.   # [cm]
        self.zr = -self.zl
        polyvals = []
        for poly in SFdata.Ez_poly:
            zil = poly.zl
            zir = poly.zr
            if zil < self.zl or zir > self.zr: 
                continue
            else:
                polyvals.append((poly.zl*1.e-2,poly.zr*1.e-2))    # all-in [m]
        polyvals = tuple(polyvals) 
        self.h = polyvals[0][1] - polyvals[0][0]
        z_parts = []
        for interval in polyvals:
            z0 = interval[0]
            z1 = z0 + self.h/4.
            z2 = z0 + self.h/2.
            z3 = z0 + (3*self.h)/4.
            z4 = z0 + self.h
            z_parts.append((z0,z1,z2,z3,z4))
        self.z_steps = NP.array(z_parts)
        return
    
    def zArray(self,n):
        return self.z_steps[n]

    def tArray(self,t0,betac):
        t1 = t0 + self.h/(4*betac)
        t2 = t0 + self.h/(2*betac)
        t3 = t0 +(3*self.h)/(4*betac)
        t4 = t0 + self.h/betac
        return NP.array([t0,t1,t2,t3,t4])

    # def z_upstream(self):
    #     return self.z_steps[0][0]

   ##   def z_downstream(self):
    #     return self.z_steps[-1][-1]

    def nsteps(self):
        return len(self.z_steps)

    def nparts(self):
        return 5
    
    def steplen(self):
        return self.h

class _DYN_G(object):
    """ DYNAC's RF-gap model """
    def __init__(self, parent):
        # _DYN_G attributes
        self.c        = PARAMS['lichtgeschwindigkeit']
        self.phis     = parent.phis
        self.freq     = parent.freq
        self.gap      = parent.gap
        self.dWf      = parent.dWf
        self.lamb     = parent.lamb
        self.EzAvg    = parent.EzAvg    #todo: check EzAvg and EzPeak again
        self.SFdata   = parent.SFdata
        self.matrix   = parent.matrix
        self.particle = parent.particle
        self.omega    = parent.omega
        self.tr       = 0.68   #todo: better use Panofski for initial value
        self.deltaW    = None
        self.particlef = None

        if self.SFdata == None:
            raise RuntimeError('_DYN_G: missing E(z) table - STOP!')
            sys.exit(1)
        # self.steps = DynacSteps(self.gap,self.SFdata)
        self.stpfac = StepFactory(self.gap,self.SFdata)
        
    def Integral1(self, z, t, h, omega, phis):
        # E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
        def E(z, t, omega=omega, phis=phis):
            z = 1.e2*z     # [cm]
            return self.SFdata.Ez0t(z,t,omega,phis)
        res = 7.*E(z[0],t[0]) + 32.*E(z[1],t[1]) +12.*E(z[2],t[2]) +32.*E(z[3],t[3]) +7.*E(z[4],t[4])
        res = res * h / 90.
        return res

    def Integral2(self, z, t, h, omega, phis):
        # E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
        def E(z, t, omega=omega, phis=phis):
            return self.SFdata.Ez0t(z,t,omega,phis)
        z = 1.e2*z     # [cm]
        res = 8.*E(z[1],t[1]) + 6.*E(z[2],t[2]) +24.*E(z[3],t[3]) + 7.*E(z[4],t[4])
        res = res * h**2 / 90.
        return res

    #todo: use beta_gamma[i], i=1,4,1
    def Integral3(self, z, t, h, bg, omega, phis):
        # E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
        def E(z, t, omega=omega, phis=phis):
            return self.SFdata.Ez0t(z,t,omega,phis)
        z = 1.e2*z     # [cm]
        bg3 = bg**3
        res = (8.*E(z[1],t[1]) + 6.*E(z[2],t[2]) + 24.*E(z[3],t[3]) + 7.*E(z[4],t[4]))/bg3
        res = res * h**2 / 90.
        return res

    def Integral4(self, z, t, h, bg, omega, phis):
        # E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
        def E(z, t, omega=omega, phis=phis):
            return self.SFdata.Ez0t(z,t,omega,phis)
        z = 1.e2*z     # [cm]
        bg3 = bg**3
        res = (2.*E(z[1],t[1]) + 3.*E(z[2],t[2]) + 18.*E(z[3],t[3]) + 7.*E(z[4],t[4]))/bg3
        res = res * h**3 / 90.
        return res

    #todo: use G1(gamma[i]), i=1,4,1
    def Jntegral1(self, z, t, h, omega, phis, gamma):
        # E = partial(self.SFdata.dEz0tdt, omega = omega, phis = phis)
        def Ep(z, t, omega=omega, phis=phis):
            return self.SFdata.dEz0tdt(z,t,omega,phis)
        def G1(gamma):
            return (gamma**2-1.)**(-1.5)/(2.*self.particle.m0c3)
        g1 = G1(gamma)
        res = (7.*Ep(z[0],t[0]) + 32.*Ep(z[1],t[1]) + 12.*Ep(z[2],t[2]) + 32.*Ep(z[3],t[3]) + 7.*Ep(z[4],t[4]))*g1
        res = res * h / 90.
        return res

    def Jntegral2(self, z, t, h, omega, phis, gamma):
        # E = partial(self.SFdata.dEz0tdt, omega = omega, phis = phis)
        def Ep(z, t, omega=omega, phis=phis):
            return self.SFdata.dEz0tdt(z,t,omega,phis)
        def G1(gamma):
            return (gamma**2-1.)**(-1.5)/(2.*self.particle.m0c3)
        g1 = G1(gamma)
        res = (8.*Ep(z[1],t[1]) + 6.*Ep(z[2],t[2]) + 24.*Ep(z[3],t[3]) + 7.*Ep(z[4],t[4]))*g1
        res = res * h**2 / 90.
        return res

    def Jntegral3(self, z, t, h, omega, phis, gamma):
        # E = partial(self.SFdata.dEz0tdt, omega = omega, phis = phis)
        def Ep(z, t, omega=omega, phis=phis):
            return self.SFdata.dEz0tdt(z,t,omega,phis)
        def G1(gamma):
            return (gamma**2-1.)**(-1.5)/(2.*self.particle.m0c3)
        g1 = G1(gamma)
        res = (2.*Ep(z[1],t[1]) + 3.*Ep(z[2],t[2]) + 18.*Ep(z[3],t[3]) + 7.*Ep(z[4],t[4]))*g1
        res = res * h**3 / 90.
        return res

    def soll_map(self, i_track):
        """ Soll mapping form (i) to (f) """
        DEBUG_SOLL = DEBUG_OFF
        particle = self.particle
        m0c2     = particle.m0c2
        beta     = particle.beta
        betac    = particle.betac
        c        = self.c
        phis     = self.phis
        omega    = self.omega
        tkin     = self.particle.tkin   # energy IN

        nsteps = self.stpfac.nsteps()    # nb-steps
        h      = self.stpfac.steplen()   # step length
    
        # z_upstream   = self.stpfac.z_upstream()
        # z_downstream = self.stpfac.z_downstream()
        # t0 = z_upstream/betac
        
        for nstep in range(nsteps):# loop steps
            # start with const beta in step
            beta  = MATH.sqrt(1.-1./(1.+tkin/m0c2)**2)
            zarr = self.stpfac.zArray(nstep)
            t0 = zarr[0]/betac
            tarr = self.stpfac.tArray(t0,betac)
            I1   = self.Integral1(zarr,tarr,h,omega,phis)
            dgamma  = I1/m0c2                               # (12)
            tkin = tkin + dgamma * m0c2
        # energy gain per rf-gap
        deltaW = tkin - particle.tkin
        # UPDATE linear NODE matrix with this deltaW
        self.matrix[EKOO, DEKOO] = deltaW
        # the parent delegates reading these properties from here
        self.tr        = deltaW/(self.EzAvg*self.gap)
        self.deltaw    = deltaW
        self.particlef = copy(particle)(tkin)     # copy is !!!IMPORTANT!!!
        # track the track
        f_track        = i_track
        f_track[Ktp.T] = tkin
        DEBUG_SOLL('SOLL',(f_track,self.tr))
        return f_track

    def map(self,i_track):
        """ Mapping from (i) to (f) """
        DEBUG_MAP = DEBUG_ON
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4]
        zp       = i_track[ZPKOO]      # [5]
        # aliases
        # ref particle
        particleS = self.particle
        m0c2      = particleS.m0c2
        m0c3      = particleS.m0c3
        # gammaS    = particleS.gamma
        # betaS     = particleS.beta
        # bgS       = particleS.gamma_beta
        # bgrootS   = MATH.sqrt(bgS)
        phiS      = self.phis
        freq      = self.freq
        omega     = self.omega
        c         = self.c
        tkinS     = particleS.tkin   # soll energy IN
        
        # particle
        converter = WConverter(tkinS,freq=freq)
        DW   = converter.Dp2pToW(zp)
        DPHI = converter.zToDphi(z)
        particle = copy(particleS)(tkinS+DW)
        gamma  = particle.gamma
        beta   = particle.beta
        betac  = particle.betac
        bg     = particle.gamma_beta
        bgroot = MATH.sqrt(bg)
        tkin   = particle.tkin
        
        # Picht transformation
        r      = NP.array([x,y])                        # (x,y)
        rp     = NP.array([xp,yp])                      # (x',y')
        R      = r*bgroot                               # (X,Y)
        Rp     = rp*bgroot+0.5*R*gamma/(gamma**2-1.)    # (X',Y')

        Dtime   = -z/(beta*c)    # time difference from ref particle IN
        # Dtime   = DPHI/omega     # time difference from ref particle IN
        t00 = self.stpfac.zArray(0)[0]/betac
        t0 = t00 + Dtime
        nsteps = self.stpfac.nsteps()        # nb-steps
        # nparts = self.steps.nparts()        # nb-parts/step
        h      = self.stpfac.steplen()       # step length
        
        for nstep in range(nsteps):         # steps loop
            # const. beta  in step
            beta  = MATH.sqrt(1.-1./(1.+tkin/m0c2)**2)
            gamma = 1./MATH.sqrt(1.- beta**2)
            bg    = beta*gamma
            K1    = omega**2/(4.*c**2*bg**3)
            betac = beta*c
            # time
            # zs,ts = self.steps(beta*c)
            # zarr = zs[nstep]
            # tarr = ts[nstep]
            zarr = self.stpfac.zArray(nstep)
            # time at z0
            # t0 = zarr[0]/betac+Dtime
            tarr = self.stpfac.tArray(t0,betac)
            I3 = self.Integral3(zarr,tarr,bg,h,omega,phiS)
            I4 = self.Integral4(zarr,tarr,bg,h,omega,phiS)
            # delta time
            Dtime = ((1.+ NP.dot(R,R)*K1)*I3 + NP.dot(R,Rp)*K1*I4)/m0c3       # (26)
            # time at z4
            t4 = t0 + Dtime + h/(betac)
            t4 = tarr[-1] + Dtime
            t0 = t4

            # energy
            I1   = self.Integral1(zarr,tarr,h,omega,phiS)
            I2   = self.Integral2(zarr,tarr,h,omega,phiS)
            # delta gamma
            Dgamma  = ((1. + NP.dot(R,R)*K1)*I1 + NP.dot(R,Rp)*K1*I2)/m0c2    # (12)
            tkin = tkin + Dgamma * m0c2
            
            # transverse
            J1   = self.Jntegral1(zarr,tarr,h,omega,phiS,gamma)
            J2   = self.Jntegral2(zarr,tarr,h,omega,phiS,gamma)
            J3   = self.Jntegral3(zarr,tarr,h,omega,phiS,gamma)
            # delta Picht
            DRp  = R*J1 + Rp*J2
            DR   = R*J2 + Rp*J3
            R    = R  + DR
            Rp   = Rp + DRp
            pass                            # end steps loop
        
        # Picht back-transformation
        particle = particle(tkin = tkin)
        gamma    = particle.gamma
        gambeta  = particle.gamma_beta
        gbroot   = MATH.sqrt(gambeta)
        rf       = R/gbroot
        rpf      = (Rp - 0.5*R*gamma/(gamma**2-1.)) / gbroot
        
        x  = rf[0]
        xp = rpf[0]
        y  = rf[1]
        yp = rpf[1]
        z  = -t4 *(betac)
        converter = WConverter(tkin,freq=freq)
        zp = converter.DWToDp2p(tkin - tkinS)

        f_track = NP.array([ x, xp, y, yp, z, zp, tkin, 1., 0., 1.])
        DEBUG_MAP('x{:+10.5f} xp{:+10.5f} y{:+10.5f} yp{:+10.5f} z {:+10.5f} zp {:+10.5f} T {:+10.5f} {:+10.5f} S{:+10.5f} {:+10.5f} '.format(f_track[0],f_track[1],f_track[2],f_track[3],f_track[4],f_track[5],f_track[6],f_track[7],f_track[8],f_track[9]))
        return f_track
            
#         
# 
# class _DYN_Gslice(object):
#     """ The DYNAC gap-model """
#     def __init__(self, parent, poly, particle):
#         # _DYN_Gslice attributes (inherited)
#         self.parent   = parent
#         self.freq     = parent.freq     # frquency
#         self.omega    = parent.omega    # Kreisfrequenz
#         self.lamb     = parent.lamb     # Wellenlaenge
#         self.SFdata   = parent.SFdata   # reference to superfish data
#         self.particle = particle        # cloned from parent
#         self.poly     = poly            # the current interval
#     #todo: calculate Tk better ??
#         # self.Tk       = 0.85
# 
#     def time_array(self,betac,h,zarr):
#         """ create arrival times azimutal positions 
#             IN-->betac = velocity
#                  h = azimutal step size
#                  zarr = azimutal positions
#             
#             OUT<--array of arrival times
#         """
#         th = h / betac         # [s] time step size (Error: was negative before)
#         # t0 is time at slice entry!
#         t0 = zarr[0] / betac   # [s] time on interval entry (Error: was negative before)
#         t4 = th+t0
#         t3 = 0.75*th+t0
#         t2 = 0.50*th+t0
#         t1 = 0.25*th+t0
#         return NP.array([t0,t1,t2,t3,t4])
#         
#     def adjust_slice_parameters(self, phin, tkin):
#         """ Adjust soll-energy-dependent parameters for this slice 
#             IN-->phin = phase of soll at entry
#                  tkin = kin.energy at entry
#         """
#         # print(phin,tkin)
#         # new soll-particle energy
#         particle   = self.particle(tkin)
#         m0c3       = particle.m0c3
#         betac      = particle.betac
# 
#         self.h = 1.e-2*2.*self.poly.dz # [m] azimutal step size
#         # z = z0 at slice entry!
#         h  = self.h
#         z0 = 1.e-2*self.poly.zl        # [m] interval left border
#         z4 = h+z0                      # [m] interval right border
#         z3 = 0.75*h+z0
#         z2 = 0.50*h+z0
#         z1 = 0.25*h+z0
#         self.zarr = NP.array([z0,z1,z2,z3,z4])                # [m]
#         tarr      = self.time_array(betac, h, self.zarr)      # [s]
#         I1        = self.Integral1(self.zarr, tarr, h, self.omega, phin, particle)
#         I3        = self.Integral3(self.zarr, tarr, h, self.omega, phin, particle)
#         # soll parameter
#         self.BETA      = particle.beta                     # beta
#         self.GAMMA     = particle.gamma                    # gamma
#         self.WIN       = tkin                              # (z0) kin.energy
#         self.PHIN      = phin                              # (z0) phase
#         self.deltaW    = I1                                # (z4) delta-kin.energy: formel (12)  E.Tanke and S.Valero
#         self.deltaT    = I3/m0c3                           # (z4) delta-time
#         self.deltaPh   = self.deltaT*self.omega            # (z4) delta-phase
#         self.WOUT      = tkin + self.deltaW                # (z4) kin.energy
#         self.PHOUT     = phin + self.deltaPh                # (z4) phase
# 
#         DEBUG_SLICE('_DYN_Gslice: {}\n'.format(self),self.__dict__)
#         DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():zarr.........[m]: ', self.zarr)
#         DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():tarr......[psec]: ', 1.e12*tarr)
#         # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():K1......[1/m**2]: ', K1)
#         DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():I1..........[MV]: ', I1)
#         # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():I2........[MV*m]: ', self.I2)
#         DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():I3........[MV*m]: ', I3)
#         # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():I4.....[MV*m**2]: ', self.I4)
#         # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():J1.........[1/m]: ', self.J1)
#         # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():J2............[]: ', self.J2)
#         # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():J3...........[m]: ', self.J3)
#         DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():tkin.......[MeV]: ', tkin)
#         DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():phin..........[]: ', phin)
#         DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():WOUT.......[MeV]: ', self.WOUT)
#         DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():PHOUT.........[]: ', self.PHOUT)
#         DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():deltaPh.......[]: ', self.deltaPh)
#         # DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():deltaT.....[sec]: ', self.deltaT)
#         DEBUG_SLICE('_DYN_Gslice:adjust_slice_parameters():deltaW.....[Kev]: ', self.deltaW*1.e3)
#         DEBUG_SLICE('============================================================ adjust_slice_parameters() end')
#         return
# 
#     def slice_map(self, i_track):
#         """ The DYNAC RF-gap model map of a slice """
#         x        = i_track[XKOO]       # [0]
#         xp       = i_track[XPKOO]      # [1]
#         y        = i_track[YKOO]       # [2]
#         yp       = i_track[YPKOO]      # [3]
#         z        = i_track[ZKOO]       # [4] z
#         zp       = i_track[ZPKOO]      # [5] dp/p
#         T        = i_track[EKOO]       # [6] summe aller dT
#         S        = i_track[SKOO]       # [8] summe aller laengen
# 
#         # particle kin.energy and phase IN
#         converter = WConverter(self.particle.T,freq=self.freq)
#         # win = (zp * (self.GAMMA+1.)/self.GAMMA +1.) * self.WIN      # falsch
#         # win = zp*(1.+1./self.GAMMA)*self.particle.T+self.WIN
#         win       = converter.Dp2pToW(zp) + self.WIN                      # kin. energy IN
#         # pin      = - z*twopi/(self.BETA*self.lamb) + self.PHIN 
#         pin       = converter.zToDphi(z)+self.PHIN                        # phase IN
#         particle  = copy(self.particle)(win)                              # particle IN
# 
#         # print('(z,zp) ({:10.5f} , {:10.5f}) win {:10.5f}'.format(i_track[Ktp.z],i_track[Ktp.zp],win))
# 
#         # aliases
#         gamma   = particle.gamma
#         beta    = particle.beta
#         gambeta = particle.gamma_beta
#         lamb    = self.lamb
#         m0c2    = particle.m0c2
#         m0c3    = particle.m0c3
#         betac   = particle.betac
#         omega   = self.omega
#         try:
#             gbroot  = MATH.sqrt(gambeta)
#         except ValueError as ex:
#             for k,v in particle.__dict__.items():
#                 print(k,v)
#             raise ex
#         gb3     = gambeta**3
#         K1      = (MATH.pi/self.lamb)**2/gb3  # [1./m**2] common factor
#  
#         # Integrale
#         tarr = self.time_array(betac, self.h, self.zarr)
#         I1   = self.Integral1(self.zarr, tarr, self.h, omega, pin, particle)
#         I2   = self.Integral2(self.zarr, tarr, self.h, omega, pin, particle)
#         I3   = self.Integral3(self.zarr, tarr, self.h, omega, pin, particle)
#         I4   = self.Integral4(self.zarr, tarr, self.h, omega, pin, particle)
#         J1   = self.Jntegral1(self.zarr, tarr, self.h, omega, pin, particle)
#         J2   = self.Jntegral2(self.zarr, tarr, self.h, omega, pin, particle)
#         J3   = self.Jntegral3(self.zarr, tarr, self.h, omega, pin, particle)
# 
#         # Picht transformation
#         r   = NP.array([x,y])                      # =(x,y)
#         rp  = NP.array([xp,yp])                    # =(x',y')
#         R   = r*gbroot                             # =(X,Y)
#         Rp  = rp*gbroot+0.5*R*gamma/(gamma**2-1.)  # =(X',Y')
#         # delta gamma
#         dgamma  = ((1. + NP.dot(R,R)*K1)*I1 + NP.dot(R,Rp)*K1*I2)/m0c2
#         deltaW  = dgamma * m0c2                    # (z4) delta-kin.energy
#         # delta time
#         dtime  = ((1. + NP.dot(R,R)*K1)*I3 + NP.dot(R,Rp)*K1*I4)/m0c3
#         dphi   = dtime * omega                     # (z4) delta-phase
#         dz1    = - beta*lamb/twopi*dphi            # z = distance from soll
#         dz2    = - betac * dtime                   # z = distance from soll
#         if abs(dz1-dz2) > 1.e-15:
#             warnings.warn('|delta-z difference| too large - should be equal')
#         DEBUG_OFF('(deltaW[KeV], dphi[mdeg]) ',(deltaW*1.e3, MATH.degrees(dphi)*1.e3))
# 
#         # mapping of reduced coordinates
#         dR     = R*J2 + Rp*J3  # delta-radius
#         dRp    = R*J1 + Rp*J2  # delta-radius primed  
#         Rf     = R + dR
#         Rpf    = Rp + dRp32.*E(z[1],t[1]) +
# 
#         # Picht transformation back
#         particle = copy(self.particle)(tkin = win + deltaW)
#         gamma    = particle.gamma
#         gambeta  = particle.gamma_beta
#         gbroot   = MATH.sqrt(gambeta)
#         rf       = Rf/gbroot
#         rpf      = (Rpf - 0.5*Rf*gamma/(gamma**2-1.)) / gbroot
# 
#         # new z
#         converter = WConverter(particle.T,freq=self.freq)
#         pout      = pin + dphi
#         dp        = pout - self.PHOUT
#         # zf        = -beta*self.lamb/twopi*dp     # z out
#         zf        = converter.DphiToz(dp)          # z OUT
#         # new dp/p 
#         wout     = win + deltaW
#         dw       = wout - self.WOUT
#         # zpf      = gamma/(gamma+1.)*dw/wout     # delta-p/p out
#         zpf      = converter.DWToDp2p(dw)         # Dp2p OUT
#   
#         # new track point
#         f_track = NP.array([x,xp,y,yp,z,zp,T,1.,S,1.])
#         f_track[XKOO]  = rf[0]
#         f_track[YKOO]  = rf[1]
#         f_track[XPKOO] = rpf[0]
#         f_track[YPKOO] = rpf[1]
#         f_track[ZKOO]  = zf
#         f_track[ZPKOO] = zpf
#         f_track[EKOO]  = T + self.deltaW
# 
#         # DEBUG_SLICE('_DYN_Gslice:slice_map():R............[m]: ', R)
#         # DEBUG_SLICE('_DYN_Gslice:slice_map():Rp.........[rad]: ', Rp)
#         # DEBUG_SLICE('_DYN_Gslice:slice_map():dR..........[um]: ', dR*1e6)
#         # DEBUG_SLICE('_DYN_Gslice:slice_map():dRp.......[urad]: ', dRp*1e6)
#         # DEBUG_SLICE('_DYN_Gslice:slice_map():dgamma..........: ', dgamma)
#         # DEBUG_SLICE('_DYN_Gslice:slice_map():dtime.....[psec]: ', dtime*1.e12)
#         # DEBUG_SLICE('_DYN_Gslice:slice_map():r=(x,y)......[m]: ', rf)
#         # DEBUG_SLICE("_DYN_Gslice:slice_map():rp=(x',y').[rad]: ", rpf)
#         # DEBUG_SLICE('============================================================ slice_map() end')
#         return f_track
# 
#     def Integral3(self, zarr, tarr, h, omega, phin, particle):
#         coeff = NP.array([0., 8., 6., 24., 7.])
#         res   = 0.
#         bg3 = particle.gamma_beta**3
#         # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
#         E = partial(self.SFdata.Ez0t, omega = omega, phi = phin)
#         for i in range(1,len(coeff)):
#             z = 1.e2*zarr[i]        # [m] --> [cm]
#     #todo: use beta_gamma[i], i=1,4,1
#             res = res + coeff[i]/bg3 * E(z, tarr[i])
#         res = res * h**2 / 90.
#         return res
# 
#     def Integral4(self, zarr, tarr, h, omega, phin, particle):
#         coeff = NP.array([0., 2., 3., 18., 7.])
#         res   = 0.
#         bg3 = particle.gamma_beta**3
#         # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
#         E = partial(self.SFdata.Ez0t, omega = omega, phi = phin)
#         for i in range(1,len(coeff)):
#             z = 1.e2*zarr[i]        # [m] --> [cm]
#     #todo: use beta_gamma[i], i=1,4,1
#             res = res + coeff[i]/bg3 * E(z, tarr[i])
#         res = res * h**3 / 90.
#         return res
# 
#     def Jntegral1(self, zarr, tarr, h, omega, phin, particle):
#         coeff = NP.array([7., 32., 12., 32., 7.])
#         res   = 0.
#         g2m1 = particle.gamma**2-1.
#     #todo: use G1(gamma[i]), i=1,4,1
#         G1    = g2m1**(-1.5)/(2.*particle.m0c3)
#         # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
#         Ep = partial(self.SFdata.dEz0tdt, omega = omega, phi = phin)
#         # print('====dE(z,0,t)/dt')
#         for i in range(0,len(coeff)):
#             z = 1.e2*zarr[i]        # [m] --> [cm]
#             res = res + coeff[i] * G1 * Ep(z, tarr[i])
#             # print('(z,Ep(z, tarr[i])) ',(z,Ep(z, tarr[i])))
#         res = res * h / 90.
#         return res
# 
#     def Jntegral2(self, zarr, tarr, h, omega, phin, particle):
#         coeff = NP.array([0., 8., 6., 24., 7.])
#         res   = 0.
#         g2m1 = particle.gamma**2-1.
#     #todo: use G1(gamma[i]), i=1,4,1
#         G1    = g2m1**(-1.5)/(2.*particle.m0c3)
#         # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
#         Ep = partial(self.SFdata.dEz0tdt, omega = omega, phi = phin)
#         for i in range(1,len(coeff)):
#             z = 1.e2*zarr[i]        # [m] --> [cm]
#             res = res + coeff[i] * G1 * Ep(z, tarr[i])
#         res = res * h**2 / 90.
#         return res
# 
#     def Jntegral3(self, zarr, tarr, h, omega, phin, particle):
#         coeff = NP.array([0., 2., 3., 18., 7.])
#         res   = 0.
#         g2m1 = particle.gamma**2-1.
#     #todo: use G1(gamma[i]), i=1,4,1
#         G1    = g2m1**(-1.5)/(2.*particle.m0c3)
#         # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
#         Ep = partial(self.SFdata.dEz0tdt, omega = omega, phi = phin)
#         for i in range(1,len(coeff)):
#             z = 1.e2*zarr[i]        # [m] --> [cm]
#             res = res + coeff[i] * G1 * Ep(z, tarr[i])
#         res = res * h**3 / 90.
#         return res
# 
def test0():
    DEBUG_TEST0 = DEBUG_ON
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
    z  = 1.0e-3
    zp = 1.0e-3
    
    i_track1 = NP.array([ 0,  0, 0,  0, 0,  0, PARAMS['sollteilchen'].tkin, 1., 0., 1.])
    i_track2 = NP.array([ x, xp, y, yp, z, zp, PARAMS['sollteilchen'].tkin, 1., 0., 1.])
    dyng = ELM.RFG(gap=0.048, SFdata=SF_tab, mapping='dyn')
    
    f_track = dyng.soll_map(i_track1)
    DEBUG_TEST0('_DYN_G.soll_map():i_track:\n', str(i_track1))
    DEBUG_TEST0('_DYN_G.soll_map():f_track:\n', str(f_track))

    DEBUG_TEST0('_DYN_G.map():i_track:\n', str(i_track2))
    tracks = []
    for i in range(10):
        f_track = dyng.map(i_track2)
        tracks.append(f_track)
        i_track2 = f_track
    DEBUG_TEST0('_DYN_G.map():f_track:')
    for track in tracks:
        print('x{:+10.5f} xp{:+10.5f} y{:+10.5f} yp{:+10.5f} z {:+10.5f} zp {:+10.5f} T {:+10.5f} {:+10.5f} S{:+10.5f} {:+10.5f} '.format(track[0],track[1],track[2],track[3],track[4],track[5],track[6],track[7],track[8],track[9]))

if __name__ == '__main__':
    test0()
