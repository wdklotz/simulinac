import sys
import numpy as np
from copy import copy
import math
from functools import partial


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

class _DYN_G(object):
    def __init__(self, parent):
        def _make_slices(parent, gap, SFdata, particle):
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

        def adjust_slice_energy(slices, phis, tkin):
            next_phase = phis
            next_tkin  = tkin
            Tklist = [0]
            for slice in slices:
                # setting phase @ slice entrance
                # slice.phis = next_phase
                # setting phase,energy @ slice entry
                slice.adjust_slice_parameters(next_phase, next_tkin)
                next_phase = slice.PHOUT    # slice OUT as next slice IN
                next_tkin  = slice.WOUT     # slice OUT as next slice IN
            deltaW  = next_tkin-tkin        # total energy kick as sum over slices
            self.tr = np.sum(np.array(Tklist))/len(Tklist)
            return deltaW
        # _DYN_G
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
                _make_slices(self, self.gap, self.SFdata, self.particle)
            DEBUG_DYN_G('_DYN_G:_make_slices()\n',list(self.slices))
            # self.deltaW = None  # initialize before using!
            # slice energy dependce
            self.deltaW = \
                adjust_slice_energy(self.slices, self.phis, self.particle.tkin)
            # update Node matrix with local deltaW
            self.matrix[EKOO, DEKOO] = self.deltaW

    def map(self,i_track):
        """Mapping from position (i) to (f) """
        f_track =  self._full_gap_map(self.slices, i_track)
        
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

    def _full_gap_map(self, slices, i_track):
        for slice in slices:
            f_track = slice._slice_map(i_track)   # map slice with DYNAC gap-model
            i_track = f_track
        return f_track

class _DYN_Gslice(object):
    """ E.Tanke, S.Valero DYNAC gap model mapping from (i) to (f) """
    def __init__(self, parent, poly, particle):
        self.parent   = parent
        self.freq     = parent.freq
#todo: set particle tkin correct at slice entry!
        self.particle = copy(particle)
        self.poly     = poly    # the current interval
        self.SFdata   = parent.SFdata
        self.lamb     = parent.lamb
        self.phis      = None  # initialized in adjust_slice_parameters
        self.WIN       = None  # initialized in adjust_slice_parameters
        self.WOUT      = None  # initialized in adjust_slice_parameters
        self.deltaW    = None  # initialized in adjust_slice_parameters
        self.PHIN      = None  # initialized in adjust_slice_parameters
        self.PHOUT     = None  # initialized in adjust_slice_parameters
        self.deltaPHI  = None  # initialized in adjust_slice_parameters

    def adjust_slice_parameters(self, phin, tkin):
        """ Adjust energy-dependent parameters for this slice """
        self.particle(tkin)
        self.phis = phin
        beta      = self.particle.beta
        lamb      = self.lamb
        
        i_track = np.array([0.,0.,0.,0.,0.,0.,0.,1.,0.,1.])
        f_track = self._slice_map(i_track)
        
        WIN   = tkin
        PHIN  = self.phis
        DW    = f_track[EKOO]
        WOUT  = WIN+DW
        DPHI  = -f_track[ZKOO]*twopi/(beta*lamb)
        PHOUT = PHIN+DPHI
        
        self.WIN      = WIN
        self.WOUT     = WOUT
        self.deltaW   = DW
        self.PHIN     = PHIN
        self.PHOUT    = PHOUT
        self.deltaPHI = DPHI
        DEBUG_SLICE('_DYN_Gslice: {}\n'.format(self),self.__dict__)
        return

    def _slice_map(self, i_track):
        def Integral1(zarr, tarr, h):
            coeff = np.array([7., 32., 12., 32., 7.])
            res   = 0.
            # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
            E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
            for i in range(0,len(coeff)):
                z = 1.e2*zarr[i]        # [m] --> [cm]
                res = res + coeff[i] * E(z, tarr[i])
            res = res * h / 90.
            return res

        def Integral2(zarr, tarr, h):
            coeff = np.array([0., 8., 6., 24., 7.])
            res   = 0.
            # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
            E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
            for i in range(1,len(coeff)):
                z = 1.e2*zarr[i]        # [m] --> [cm]
                res = res + coeff[i] * E(z, tarr[i])
            res = res * h**2 / 90.
            return res
#todo: use beta_gamma[i], i=1,4,1
        def Integral3(zarr, tarr, h):
            coeff = np.array([0., 8., 6., 24., 7.])
            res   = 0.
            # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
            E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
            for i in range(1,len(coeff)):
                z = 1.e2*zarr[i]        # [m] --> [cm]
                res = res + coeff[i]/beta_gamma**3 * E(z, tarr[i])
            res = res * h**2 / 90.
            return res

#todo: use beta_gamma[i], i=1,4,1
        def Integral4(zarr, tarr, h):
            coeff = np.array([0., 2., 3., 18., 7.])
            res   = 0.
            # self.SFdata:Ez0t(self,z,t,omega,phi): time dependent field value at location z
            E = partial(self.SFdata.Ez0t, omega = omega, phi = phis)
            for i in range(1,len(coeff)):
                z = 1.e2*zarr[i]        # [m] --> [cm]
                res = res + coeff[i]/beta_gamma**3 * E(z, tarr[i])
            res = res * h**3 / 90.
            return res

#todo: use G1(gamma[i]), i=1,4,1
        def Jntegral1(zarr, tarr, h):
            coeff = np.array([7., 32., 12., 32., 7.])
            res   = 0.
            G1    = g2m1**(-1.5)/2./m0c3
            # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
            Ep = partial(self.SFdata.dEz0tdt, omega = omega, phi = phis)
            for i in range(0,len(coeff)):
                z = 1.e2*zarr[i]        # [m] --> [cm]
                res = res + coeff[i] * G1 * Ep(z, tarr[i])
            res = res * h / 90.
            return res

#todo: use G1(gamma[i]), i=1,4,1
        def Jntegral2(zarr, tarr, h):
            coeff = np.array([0., 8., 6., 24., 7.])
            res   = 0.
            G1    = g2m1**(-1.5)/2./m0c3
            # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
            Ep = partial(self.SFdata.dEz0tdt, omega = omega, phi = phis)
            for i in range(1,len(coeff)):
                z = 1.e2*zarr[i]        # [m] --> [cm]
                res = res + coeff[i] * G1 * Ep(z, tarr[i])
            res = res * h**2 / 90.
            return res

#todo: use G1(gamma[i]), i=1,4,1
        def Jntegral3(zarr, tarr, h):
            coeff = np.array([0., 2., 3., 18., 7.])
            res   = 0.
            G1    = g2m1**(-1.5)/2./m0c3
            # self.SFdata.dEz0tdt(self,z,t,omega,phi): time derivative of field value at location z"""
            Ep = partial(self.SFdata.dEz0tdt, omega = omega, phi = phis)
            for i in range(1,len(coeff)):
                z = 1.e2*zarr[i]        # [m] --> [cm]
                res = res + coeff[i] * G1 * Ep(z, tarr[i])
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
        c           = PARAMS['lichtgeschwindigkeit']
        beta        = self.particle.beta
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
        K1          = (math.pi/lamb)**2/beta_gamma**3  # [1./m**2]

        # Picht transformation
        r   = np.array([x,y])                      # =(x,y)
        rp  = np.array([xp,yp])                    # =(x',y')
        R   = r*gbroot                             # =(X,Y)
        Rp  = rp*gbroot+0.5*R*gamma/(gamma**2-1.)  # =(X',Y')

        h  = 1.e-2*2.*self.poly.dz    # [m] azimutal step size
        th = h/betac                  # [s] time step size
        # z at slice entry!
        z0 = 1.e-2*self.poly.zl       # [m] interval left border
        z4 = h+z0                     # [m] interval right border
        z3 = 0.75*h+z0
        z2 = 0.50*h+z0
        z1 = 0.25*h+z0
        # t0 at slice entry!
        t0 = - z0 / betac             # [s]
        t4 = th+t0
        t3 = 0.75*th+t0
        t2 = 0.50*th+t0
        t1 = 0.26*th+t0
        zarr = np.array([z0,z1,z2,z3,z4])   # [m]
        tarr = np.array([t0,t1,t2,t3,t4])   # [s]

        I1      = Integral1(zarr, tarr, h)
        I2      = Integral2(zarr, tarr, h)
        dgamma  = ((1. + np.dot(R,R)*K1)*I1 + np.dot(R,Rp)*K1*I2)/m0c2
        deltaW  = dgamma * m0c2                                   # delta kin. energy
#todo: normalize dp/p to W @ (f)?
        dp2p    = gamma/(gamma+1)*deltaW/W                        # dp/p

        I3     = Integral3(zarr, tarr, h)
        I4     = Integral4(zarr, tarr, h)
        dtime  = ((1. + np.dot(R,R)*K1)*I3 + np.dot(R,Rp)*K1*I4)/m0c3
        dphi   = omega * dtime                                     # dphi = phase jump
        dz0    = - beta*lamb/twopi*dphi                            # z = distance from soll

        J1     = Jntegral1(zarr, tarr, h)
        J2     = Jntegral2(zarr, tarr, h)
        J3     = Jntegral3(zarr, tarr, h)
        dR     = R*J2 + Rp*J3
        dRp    = R*J1 + Rp*J2
        Rf     = R + dR        # reduzierte Koordinaten @ (f)
        Rpf    = Rp + dRp      # reduzierte Koordinaten @ (f)
#todo: use beta*gamma (f)
        # Picht back-transformation
        xyf  = Rf/gbroot
        xypf = (Rpf - 0.5*Rf*gamma/(gamma**2-1.)) / gbroot

        f_track = np.array([x,xp,y,yp,z,zp,T,1.,s,1.])
        f_track[XKOO]  = xyf[0]
        f_track[YKOO]  = xyf[1]
        f_track[XPKOO] = xypf[0]
        f_track[YPKOO] = xypf[1]
        f_track[ZKOO]  = dz0
        f_track[ZPKOO] = dp2p
        f_track[EKOO]  = T + deltaW

        DEBUG_SLICE('_DYN_Gslice:_slice_map():zarr.........[m]: ', zarr)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():tarr......[psec]: ', 1.e12*tarr)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():K1......[1/m**2]: ',K1)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():R............[m]: ', R)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():Rp.........[rad]: ', Rp)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():dR..........[um]: ', dR*1e6)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():dRp.......[urad]: ', dRp*1e6)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():I1..........[MV]: ', I1)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():I2........[MV*m]: ', I2)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():I3........[MV*m]: ', I3)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():I4.....[MV*m**2]: ', I4)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():J1.........[1/m]: ', J1)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():J2............[]: ', J2)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():J3...........[m]: ', J3)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():dgamma..........: ', dgamma)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():deltaW.....[Kev]: ', deltaW*1.e3)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():dtime.....[psec]: ', dtime*1.e12)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():dphase....[mdeg]: ', math.degrees(dphi)*1.e3)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():z ..........[um]: ', dz0*1e6)
        DEBUG_SLICE('_DYN_Gslice:_slice_map():r=(x,y)......[m]: ', xyf)
        DEBUG_SLICE("_DYN_Gslice:_slice_map():rp=(x',y').[rad]: ", xypf)
        DEBUG_SLICE('============================================================ slice end')

        return f_track

def test0():
    import elements as ELM
    # from tracks import Track

    print('-----------------------------------TEST 0----------------')
    print('test _DYN_Gslice:_slice_map()...')
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
    DEBUG_TEST0('_DYN_Gslice:test0():i_track:\n', str(i_track))
    f_track = dyng.map(i_track)
    DEBUG_TEST0('_DYN_Gslice:test0():f_track:\n', str(f_track))

if __name__ == '__main__':
    test0()
