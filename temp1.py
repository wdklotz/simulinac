import sys
import numpy as np
from copy import copy
import math
from functools import partial


from setutil import DEBUG, arrprnt, PARAMS
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO
from Ez0 import SFdata, Ipoly

def DEBUG_ON(string,arg = '',end = '\n'):
    DEBUG(string,arg,end)
def DEBUG_OFF(string,arg = '',end = '\n'):
    pass
# DEBUG__*
DEBUG_MAP    = DEBUG_ON

twopi = 2.*math.pi

class _DYN_G(object):
    def __init__(self, parent):
        def _make_slices():
            """Slice the RF gap"""
            slices = []
            zl = -self.gap/2.*100.   # [m] --> [cm]
            zr = -zl
            E0z = 0.
            z = 0.
            for poly in self.SFdata.Ez_poly:
                zil = poly.zl
                zir = poly.zr
                if zil < zl or zir > zr: continue
                slice = _DYN_Gslice(self, poly, self.particle)  # instanciate _DYN_Gslices
                slices.append(slice)
            return slices
        
        def adjust_slice_energy():
            pass

        self.particle = parent.particle
        self.gap      = parent.gap
        self.lamb     = parent.lamb
        self.freq     = parent.freq
        self.phis     = parent.phis
        self.SFdata   = parent.SFdata
        self.deltaW   = None  # initialize before using!

        if self.SFdata == None:
            raise RuntimeError('_DYN_G: missing E(z) table - STOP!')
            sys.exit(1)
        else:
            self.slices = _make_slices()
            DEBUG_MAP('_DYN_G:_make_slices()\n',list(self.slices))
            # self.deltaW = adjust_slice_energy()
            # self.matrix[EKOO, DEKOO] = self.deltaW
            # if DEBUG_MAP == DEBUG_ON: print(self.slices)
    
    def _full_gap_map(self,i_track):
        for slice in self.slices[1:2]:
            f_track = slice._slice_map(i_track)
        return f_track

    def map(self,i_track):
        return self._full_gap_map(i_track)
        
class _DYN_Gslice(object):
    """
    E.Tanke, S.Valero DYNAC gap model mapping from (i) to (f)
    """
    def __init__(self, parent, poly, particle):
        self.parent   = parent
        self.particle = copy(particle)
        self.poly     = poly    # the current interval
        self.SFdata   = parent.SFdata
        self.lamb     = parent.lamb
        self.freq     = parent.freq
        self.phis     = parent.phis
        
    def _slice_map(self, i_track):
        def Integral1(zarr, tarr, h):
            coeff = np.array([7., 32., 12., 32., 7.])
            res   = 0.
            # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
            E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
            for i in range(0,len(coeff)):
                res = res + coeff[i] * E(zarr[i], tarr[i])
            res = res * h / 90.
            return res

        def Integral2(zarr, tarr, h):
            coeff = np.array([0., 8., 6., 24., 7.])
            res   = 0.
            # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
            E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
            for i in range(1,len(coeff)):
                res = res + coeff[i] * E(zarr[i], tarr[i])
            res = res * h**2 / 90.
            return res
#todo: use beta_gamma[i], i=1,4,1
        def Integral3(zarr, tarr, h):
            coeff = np.array([0., 8., 6., 24., 7.])
            res   = 0.
            # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
            E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
            for i in range(1,len(coeff)):
                res = res + coeff[i]/beta_gamma**3 * E(zarr[i], tarr[i])
            res = res * h**2 / 90.
            return res

#todo: use beta_gamma[i], i=1,4,1
        def Integral4(zarr, tarr, h):
            coeff = np.array([0., 2., 3., 18., 7.])
            res   = 0.
            # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
            E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
            for i in range(1,len(coeff)):
                res = res + coeff[i]/beta_gamma**3 * E(zarr[i], tarr[i])
            res = res * h**3 / 90.
            return res

#todo: use G1(gamma[i]), i=1,4,1
        def Jntegral1(zarr, tarr, h):
            coeff = np.array([7., 32., 12., 32., 7.])
            res   = 0.
            G1    = g2m1**(-1.5)/2./m0c3
            # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
            E = partial(self.SFdata.dEz0tdt, omega = omega, phi = phis)
            for i in range(0,len(coeff)):
                res = res + coeff[i] * G1 * E(zarr[i], tarr[i])
            res = res * h / 90.
            return res

#todo: use G1(gamma[i]), i=1,4,1
        def Jntegral2(zarr, tarr, h):
            coeff = np.array([0., 8., 6., 24., 7.])
            res   = 0.
            G1    = g2m1**(-1.5)/2./m0c3
            # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
            E = partial(self.SFdata.dEz0tdt, omega = omega, phi = phis)
            for i in range(1,len(coeff)):
                res = res + coeff[i] * G1 * E(zarr[i], tarr[i])
            res = res * h**2 / 90.
            return res

#todo: use G1(gamma[i]), i=1,4,1
        def Jntegral3(zarr, tarr, h):
            coeff = np.array([0., 2., 3., 18., 7.])
            res   = 0.
            G1    = g2m1**(-1.5)/2./m0c3
            # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
            E = partial(self.SFdata.dEz0tdt, omega = omega, phi = phis)
            for i in range(1,len(coeff)):
                res = res + coeff[i] * G1 * E(zarr[i], tarr[i])
            res = res * h**3 / 90.
            return res

        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z
        zp       = i_track[ZPKOO]      # [5] dp/p
        T        = i_track[EKOO]       # [6] summe aller delta-T
        s        = i_track[SKOO]       # [8] summe aller laengen

        W           = self.particle.tkin
        beta        = self.particle.beta
        c           = PARAMS['lichtgeschwindigkeit']
        betac       = beta*c
        gamma       = self.particle.gamma
        beta_gamma  = self.particle.gamma_beta
        m0c2        = self.particle.e0
        m0c3        = m0c2*c
        omega       = twopi*self.freq
        phis        = self.phis
        lamb        = self.lamb
        g2m1        = gamma**2-1.
        gbroot      = math.sqrt(beta_gamma)
        # Picht transformation
        r   = np.array([x,y])                      # =(x,y)
        rp  = np.array([xp,yp])                    # =(x',y')
        R   = r*gbroot                             # =(X,Y)
        Rp  = rp*gbroot+0.5*R*gamma/(gamma**2-1.)  # =(X',Y')
        K1  = (math.pi/lamb)**2/beta_gamma**3      # [1./m**2]

        h  = 1.e-2*2.*self.poly.dz    # [m] azimutal step size
        th = h/betac                  # [s] time step size
        z0 = 1.e-2*self.poly.zl       # [m]
        z4 = h+z0
        z3 = 0.75*h+z0
        z2 = 0.50*h+z0
        z1 = 0.25*h+z0
        t0 = 0.
        t4 = th+t0
        t3 = 0.75*th+t0
        t2 = 0.50*th+t0
        t1 = 0.26*th+t0
        zarr = np.array([z0,z1,z2,z3,z4])
        tarr = np.array([t0,t1,t2,t3,t4])

        intg1   = Integral1(zarr, tarr, h)
        intg2   = Integral2(zarr, tarr, h)
        dgamma  = ((1. + np.dot(R,R)*K1)*intg1 + np.dot(R,Rp)*K1*intg2)/m0c2
        deltaW  = dgamma * m0c2                                   # delta kin. energy
        dp2p    = gamma/(gamma+1)*deltaW/W                        # dp/p
        self.parent.deltaW = deltaW                               # parent property

        intg3  = Integral3(zarr, tarr, h)
        intg4  = Integral4(zarr, tarr, h)
        dtime  = ((1. + np.dot(R,R)*K1)*intg3 + np.dot(R,Rp)*K1*intg4)/m0c3    # time difference
        dphi   = omega * dtime                                     # dphi = phase jump
        dz0    = - beta*lamb/twopi*dphi                            # z = distance from soll

        jntg1  = Jntegral1(zarr, tarr, h)
        jntg2  = Jntegral2(zarr, tarr, h)
        jntg3  = Jntegral3(zarr, tarr, h)
        dR     = R*jntg2 + Rp*jntg3
        dRp    = R*jntg1 + Rp*jntg2
        Rf     = R + dR
        Rpf    = Rp + dRp

        f_track = i_track
        f_track[ZKOO]  = f_track[ZKOO] + dz0
        f_track[ZPKOO] = dp2p
        f_track[EKOO]  = f_track[EKOO] + deltaW

        DEBUG_MAP('_DYN_Gslice:_slice_map():zarr.....[m]:\n', list(zarr))
        DEBUG_MAP('_DYN_Gslice:_slice_map():tarr..[psec]:\n', list(map(lambda x: x*1.e12,tarr)))
        DEBUG_MAP('_DYN_Gslice:_slice_map():K1..[1/m**2]: ',K1)
        DEBUG_MAP('_DYN_Gslice:_slice_map():R........[m]: ', R)
        DEBUG_MAP('_DYN_Gslice:_slice_map():Rp.....[rad]: ', Rp)
        DEBUG_MAP('_DYN_Gslice:_slice_map():dR......[um]: ', dR*1e6)
        DEBUG_MAP('_DYN_Gslice:_slice_map():dRp...[urad]: ', dRp*1e6)
        DEBUG_MAP('_DYN_Gslice:_slice_map():I1......[MV]: ', intg1)
        DEBUG_MAP('_DYN_Gslice:_slice_map():I2....[MV*m]: ', intg2)
        DEBUG_MAP('_DYN_Gslice:_slice_map():I3....[MV*m]: ', intg3)
        DEBUG_MAP('_DYN_Gslice:_slice_map():I4 [MV*m**2]: ', intg4)
        DEBUG_MAP('_DYN_Gslice:_slice_map():J1.....[1/m]: ', jntg1)
        DEBUG_MAP('_DYN_Gslice:_slice_map():J2........[]: ', jntg2)
        DEBUG_MAP('_DYN_Gslice:_slice_map():J3.......[m]: ', jntg3)
        DEBUG_MAP('_DYN_Gslice:_slice_map():dgamma......: ', dgamma)
        DEBUG_MAP('_DYN_Gslice:_slice_map():deltaW [Kev]: ', deltaW*1.e3)
        DEBUG_MAP('_DYN_Gslice:_slice_map():dtime [psec]: ', dtime*1.e12)
        DEBUG_MAP('_DYN_Gslice:_slice_map():dphase [deg]: ', math.degrees(dphi))
        DEBUG_MAP('_DYN_Gslice:_slice_map():z ......[um]: ', dz0*1e6)

        return f_track

    # def map(self, i_track):
    #     """Mapping from position (i) to (f)"""
    #     f_track = self._mapping(i_track)   # NOTE: use local map with sliced TTFGap
    #     self.particlef = copy(self.particle)(f_track[EKOO])

    # for DEBUGGING
    #     if DEBUG_MAP == DEBUG_ON:
    #         f = f_track.copy()
    #         for i in range(len(f_track)-4):
    #             f[i] =f[i]*1.e3
    #         arrprnt(f,fmt='{:6.3g},',txt='ttf_map: ')
    #     return f_track
def test0():
    import elements as ELM
    # from tracks import Track
    
    print('-----------------------------------TEST 0----------------')
    input_file='SF_WDK2g44.TBL'
    Epeak = PARAMS['Ez_feld']*1.8055 # [Mv/m] Epeak/Eav fuer INTG(NG(von 0 bis 2.2*sigma)
    SF_tab = SFdata(input_file,Epeak)
    x  = 1.0e-2
    xp = 0.0e-3
    y  = 0.0e-2
    yp = 2.0e-3
    z  = 0.0
    zp = 0.0
    i_track = np.array([ x, xp, y, yp, z, zp, PARAMS['sollteilchen'].tkin, 1., 0., 1.])

    dyng = ELM.RFG(gap = 0.048,SFdata = SF_tab,mapping = 'dyn')
    DEBUG_MAP('_DYN_Gslice:test0():i_track:\n', str(i_track))
    f_track = dyng.map(i_track)
    DEBUG_MAP('_DYN_Gslice:test0():f_track:\n', str(f_track))

if __name__ == '__main__':
    test0()
