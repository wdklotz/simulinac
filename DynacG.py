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
import numpy as NP
from copy import copy
import math as MATH
from functools import partial
import warnings

from setutil import DEBUG, arrprnt, PARAMS, tblprnt, Ktp, WConverter
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO
from Ez0 import SFdata, Ipoly

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
        self.stpfac = StepFactory(self.gap,self.SFdata)
        
    def Integral1(self, z, t, h, omega, phis):
        #todo: use non_const beta, i.e. gamma[i], i=1,4,1 in inegrals
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
        #todo: don't use do_step(), make faster!
        """ Soll mapping form (i) to (f) """
        DEBUG_SOLL = DEBUG_OFF
        tkin     = self.particle.tkin
        m0c2     = PARAMS['proton_mass']
        phiS     = self.phis
        omega    = self.omega

        tkIN = tkin   # energy IN
        nsteps = self.stpfac.nsteps()    # nb-steps
    
        for nstep in range(nsteps):# loop steps
            # do with const beta in step
            gamma = 1.+ tkin/m0c2
            R = Rp = NP.array([0,0])
            z = 0
            DR,DRp,Dgamma,Dtime = self.do_step(nstep,z,gamma,R,Rp,omega,phiS)
            tkin = tkin + Dgamma*m0c2
        # energy gain per rf-gap
        deltaW = tkin - tkIN
        # UPDATE linear NODE matrix with this deltaW
        self.matrix[EKOO, DEKOO] = deltaW
        # the parent delegates reading these properties from here
        self.tr        = deltaW/(self.EzAvg*self.gap)
        self.deltaw    = deltaW
        self.particlef = copy(self.particle)(tkin)  # copy is !!!IMPORTANT!!!
        # track the track
        f_track = i_track
        f_track[Ktp.T] += deltaW
        DEBUG_SOLL('SOLL',(f_track,self.tr))
        return f_track

    def do_step(self,nstep,z,gamma,R,Rp,omega,phiS):
        """
        Full step through the DYNAC intervall with its 4 parts
        """
        c = PARAMS['lichtgeschwindigkeit']
        h = self.stpfac.steplen()  
        m0c2 = PARAMS['proton_mass']
        m0c3 = m0c2*c     
        bg = MATH.sqrt(gamma**2-1)
        beta = bg/gamma
        betac = beta*c

        zarr = self.stpfac.zArray(nstep)
        t0   = (zarr[0]-z)/betac
        tarr = self.stpfac.tArray(t0,betac)

        I1   = self.Integral1(zarr,tarr,h,omega,phiS)
        I2   = self.Integral2(zarr,tarr,h,omega,phiS)
        I3   = self.Integral3(zarr,tarr,h,bg,omega,phiS)
        I4   = self.Integral4(zarr,tarr,h,bg,omega,phiS)
        J1   = self.Jntegral1(zarr,tarr,h,omega,phiS,gamma)
        J2   = self.Jntegral2(zarr,tarr,h,omega,phiS,gamma)
        J3   = self.Jntegral3(zarr,tarr,h,omega,phiS,gamma)
        K1   = omega**2/(4.*c**2*bg**3)
        Dgamma  = ((1.+ NP.dot(R,R)*K1)*I1 + NP.dot(R,Rp)*K1*I2)/m0c2     # (12)
        Dtime   = ((1.+ NP.dot(R,R)*K1)*I3 + NP.dot(R,Rp)*K1*I4)/m0c3     # (26)
        DRp     = R*J1 + Rp*J2                                            # (33)
        DR      = R*J2 + Rp*J3                                            # (38)
        return DR,DRp,Dgamma,Dtime

    def map(self,i_track):
        """ Mapping from (i) to (f) """
        DEBUG_MAP = DEBUG_OFF
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4]
        zp       = i_track[ZPKOO]      # [5]
        S        = i_track[SKOO]

        # aliases
        c      = PARAMS['lichtgeschwindigkeit']
        h      = self.stpfac.steplen()  
        m0c2   = PARAMS['proton_mass']
        phiS   = self.phis
        freq   = self.freq
        omega  = self.omega
        nsteps = self.stpfac.nsteps()        # nb-steps

        # ref particle
        tkinS     = self.particle.tkin   # soll energy IN
        gammaS    = 1.+ tkinS/m0c2
        bgS       = MATH.sqrt(gammaS**2-1)
        betaS     = bgS/gammaS
        zS        = 0
       
        # particle
        converter = WConverter(tkinS,freq=freq)
        DW        = converter.Dp2pToW(zp)
        tkin      = tkinS+DW
        gamma     = 1.+ tkin/m0c2
        bg        = MATH.sqrt(gamma**2-1)
        
        # Picht transformation
        # ref. particle
        RS = RpS = NP.array([0,0])
        # particle
        bgroot = MATH.sqrt(bg)
        r      = NP.array([x,y])                        # (x,y)
        rp     = NP.array([xp,yp])                      # (x',y')
        R      = r*bgroot                               # (X,Y)
        Rp     = rp*bgroot+0.5*R*gamma/(gamma**2-1.)    # (X',Y')

        for nstep in range(nsteps):         # steps loop
            # ref particle step
            #todo: replace by sigle liner
            gammaS = 1.+ tkinS/m0c2
            DRS,DRpS,DgammaS,DtimeS = self.do_step(nstep,zS,gammaS,RS,RpS,omega,phiS)
            bgS = MATH.sqrt(gammaS**2-1)
            betaS = bgS/gammaS
            # ref. particle time at t4
            DtimeS = h/(betaS*c)
            # ref. particle energy
            tkinS = tkinS + DgammaS*m0c2
            gammaS = 1.+ tkinS/m0c2

            # particle step
            DR,DRp,Dgamma,Dtime = self.do_step(nstep,z,gamma,R,Rp,omega,phiS)
            bg    = MATH.sqrt(gamma**2-1)
            beta  = bg/gamma
            # partile time at z4
            DtimeP = h/(beta*c) + Dtime
            # correction to z from this step   (der Knackpunt der mich 1 Woche Arbeit kostete!)
            z = z + (DtimeP - DtimeS)*(beta*c)
            # particle energy
            tkin = tkin + Dgamma*m0c2
            gamma = 1.+ tkin/m0c2

            # transverse
            R    = R  + DR
            Rp   = Rp + DRp
            pass                            # end steps loop
        
        # Picht back-transformation
        gamma    = 1.+ tkin/m0c2
        bg       = MATH.sqrt(gamma**2-1)
        gbroot   = MATH.sqrt(bg)
        rf       = R/gbroot
        rpf      = (Rp - 0.5*R*gamma/(gamma**2-1.)) / gbroot
        
        x  = rf[0]
        xp = rpf[0]
        y  = rf[1]
        yp = rpf[1]
        converter = WConverter(tkinS,freq=freq)
        zp = converter.DWToDp2p(tkin-tkinS)

        f_track = NP.array([ x, xp, y, yp, z, zp, tkin, 1., S, 1.])
        # DEBUG_MAP('MAP',(f_track,None))
        DEBUG_MAP('MAP ','x{:+10.5f} xp{:+10.5f} y{:+10.5f} yp{:+10.5f} z {:+10.5f} zp {:+10.5f} T {:+10.5f} {:+10.5f} S{:+10.5f} {:+10.5f} '.format(f_track[0],f_track[1],f_track[2],f_track[3],f_track[4],f_track[5],f_track[6],f_track[7],f_track[8],f_track[9]))
        return f_track
            
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
