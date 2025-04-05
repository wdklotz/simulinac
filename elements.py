__version__='11.0.2.4'
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
import unittest
from separatrix import w2phi
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, DSKOO, MDIM
from setutil import DEBUG_ON,DEBUG_OFF,Proton,PARAMS,wrapRED
from setutil import Ktp
import warnings
import math    as M
import numpy   as NP
import setutil as UTIL

twopi = 2.*M.pi
NP.set_printoptions(linewidth = 132, formatter = {'float': '{:>8.5g}'.format}) # numpy pretty printing
def K(gradient, particle):
    """ quad strength K[1/m**2] for protons, gradient[T/m] """
    return 0.31952 * gradient/particle.gamma_beta

""" ------- The mother of all lattice element objects (a.k.a. nodes)# ------ """
class Node(object):
    """ base class for transfer matrices (linear map)
        ii)  is a dictionary (DictObject base class)
        ii)  each instance holds its copy of the refrence particle (self.particle)
    """
    def __init__(self, tsoll):
        self.type         = self.__class__.__name__  # self's node type
        self._tsoll       = tsoll # private read-only
        self.accelerating = False
        self.label        = ''
        self.viseo        = 0                   # default - invisible
        self.particle     = Proton(self._tsoll) # Sollteilchen
        self.particlef    = Proton(self._tsoll)                 # Teichen am Ausgang
        self.matrix       = NP.eye(MDIM)        # MDIMxMDIM zero matrix used here
        self.position     = (0,0,0)             # [entrance, middle, exit]
        self.length       = 0.
        self.aperture     = None
        self.next         = self      # next Node
        self.prev         = self      # pevious Node
        self.twiss        = None      # twiss functions @ middle of Node
        self.sigxy        = None      # envelope function @ middle of Node
        self.sec          = ''        # section
        self.soll_track   = None      # soll track @ exit of Node

    @property
    def tsoll(self): return self._tsoll    # Sollenergie

    def submx(self, n = MDIM, m = MDIM):
        # return upper left n, m submatrix
        return self.matrix[:n, :m]
    def map(self,i_track):
        """ standard mapping with T3D matrix """
        f_track = NP.dot(self.matrix,i_track)
        return f_track
    def toString(self):
        s = repr(self)
        for k,v in self.__dict__.items():
            s+= '\n{}:{}'.format(k,v)
        return s
    def adjust_energy(self, tkin):
        """ minimal adjust """
        self._tsoll = tkin
        self.particle(tkin)
        self.particlef(tkin)
        return self
    def shorten(self, length):
        """ nothing to shorten """
        return self
    def make_slices(self, anz=1):
        """ nothing to slice """
        slices = [self]
        # if self.length == 0. or anz <=0:    anz = 1
        # step = self.length/anz 
        # mx   = self.shorten(step)
        # for i in range(anz):
        #     slices.append(mx)
        return slices
    def aper_check(self,new_tp,s,**kwargs):
        """ nothing to check """
        return False   

class I(Node):
    """  Unit matrix Node """
    def __init__(self,label,tsoll):
        super().__init__(tsoll)
        self.label    = label
class MRK(I):
    """ 
    Marker node (a.k.a element): Each marker is parent of an agent that does the specific action.
    The action can be bypassed if the 'maction'-FLAG is False.
    """
    def __init__(self, label, active, viseo, tsoll):
        super().__init__(label,tsoll)
        self.active     = active
        self.viseo      = viseo if self.active else 0.
    def no_action(self,*args): pass
class D(Node):
    """  Trace3D drift space  """
    def __init__(self, label, length, aperture, tsoll):
        super().__init__(tsoll)
        self.label    = label
        self.length   = length
        self.aperture = aperture
        self.viseo    = 0.
        g = self.particle.gamma
        self.matrix[XKOO, XPKOO] = self.matrix[YKOO, YPKOO] = self.length
        self.matrix[ZKOO, ZPKOO] = self.length/(g*g)
        self.matrix[SKOO, DSKOO] = self.length # Delta-s longitudinal length increase

    def adjust_energy(self, tkin):
        adjusted = D(self.label, self.length, self.aperture, tkin)
        return adjusted
    def shorten(self, length):
        shortend =  D(self.label, length, self.aperture, self.tsoll)
        return shortend
class DKD(D):
    """  Trace3D drift spacer for Drift-Kick-Drift sandwich (for use with DYNAC CAVNUM cavities) """
    def __init__(self, label, length, aperture, tsoll):
        super().__init__(label, length, aperture, tsoll)

    def adjust_energy(self, tkin):
        adjusted = DKD(self.label, self.length, self.aperture, tkin)
        return adjusted
    def shorten(self, length):
        shortend =  DKD(self.label, length, self.aperture, self.tsoll)
        return shortend
class QF(Node):
    """ Trace3D focussing quad """
    def __init__(self, label, grad, length, aperture, tsoll):
        super().__init__(tsoll)
        self.label    = label
        self.grad     = abs(grad)                        # [T/m]
        self.length   = length
        self.aperture = aperture
        self.viseo    = +0.5
        self.thin     = False
        self.matrix   = self.mx()
    def mx(self):
        m = NP.eye(MDIM)
        g = self.particle.gamma
        rzz12 = self.length/(g*g)
        k02 = K(self.grad,self.particle)       # [1/m**2]
        k0w   = M.sqrt(k02)
        l = self.length
        phi = l*k0w
        f = 1./(k02*l)
        if self.thin != True:
            """ thick quad """
            # focusing
            cf   = M.cos(phi)
            sf   = M.sin(phi)/k0w
            cfp  = -k0w*M.sin(phi)
            sfp  = cf 
            # defocusing
            cd  =  M.cosh(phi)
            sd  =  M.sinh(phi)/k0w
            cdp =  k0w*M.sinh(phi)
            sdp = cd
            if self.type == "QF":
                m[XKOO, XKOO]  = cf; m[XKOO, XPKOO] = sf; m[XPKOO, XKOO] = cfp; m[XPKOO, XPKOO] = sfp
                m[YKOO, YKOO]  = cd; m[YKOO, YPKOO] = sd; m[YPKOO, YKOO] = cdp; m[YPKOO, YPKOO] = sdp
            else:
                m[XKOO, XKOO]  = cd; m[XKOO, XPKOO] = sd; m[XPKOO, XKOO] = cdp; m[XPKOO, XPKOO] = sdp
                m[YKOO, YKOO]  = cf; m[YKOO, YPKOO] = sf; m[YPKOO, YKOO] = cfp; m[YPKOO, YPKOO] = sfp
        else:
            """ thin quad D-Kick-D"""
            if self.type == "QF":
                pass
            else:
                f = -f
            m[XKOO, XKOO]  = 1.-l/(2.*f); m[XKOO, XPKOO] = l-l**2/(6*f); m[XPKOO, XKOO] = -1./f; m[XPKOO, XPKOO] = m[XKOO, XKOO]
            f = -f
            m[YKOO, YKOO]  = 1.-l/(2.*f); m[YKOO, YPKOO] = l-l**2/(6*f); m[YPKOO, YKOO] = -1./f; m[YPKOO, YPKOO] = m[YKOO, YKOO]
        m[ZKOO, ZPKOO] = rzz12
        m[SKOO, DSKOO]  = self.length # length increase
        return m
    def adjust_energy(self, tkin):
        adjusted = QF(self.label, self.grad, self.length, self.aperture, tkin)
        return adjusted
    def shorten(self, length):
        shortened = QF(self.label, self.grad, length, self.aperture, self.tsoll)
        return shortened
    def aper_check(self,new_tp,s,**kwargs):
        new_point=new_tp()
        fifo_xy= kwargs['fifo_xy']
        sfifo_xy=kwargs['sfifo_xy']
        lost=False
        # transverse apertures
        if self.aperture != None and not (abs(new_point[Ktp.x]) < self.aperture or abs(new_point[Ktp.y]) < self.aperture):
            fifo_xy.append(f'loss (x|y) ({new_point[Ktp.x]:.3e},{new_point[Ktp.y]:.3e}) at {s:.4e} m')
            sfifo_xy.append(s)
            lost = True
        return lost
class QD(QF):
    """ 
    Trace3D defocussing quad
    """
    def __init__(self, label, grad, length, aperture, tsoll):
        super().__init__(label, grad, length, aperture, tsoll)
        self.viseo = -0.5
    def adjust_energy(self, tkin):
        adjusted = QD(self.label, self.grad, self.length, self.aperture, tkin)
        return adjusted
    def shorten(self, length):
        shortened = QD(self.label, self.grad, length,  self.aperture, self.tsoll)
        return shortened
class RFG(Node):   
    """  RF-gap of zero length for different kick gap-models """
    def __init__(self,label, tsoll):
        super().__init__(tsoll)
        self.accelerating = True
        self.cavlen       = None
        self.deltaW       = None
        self.dWf          = UTIL.FLAGS['dWf']
        self.EzPeak       = None
        self.freq         = None
        self.gap          = None
        self.label        = label
        self.lamb         = None
        self.mapper       = None
        self.mapping      = UTIL.FLAGS['mapping']
        self.omega        = None
        self.phisoll      = None
        self.SFdata       = None
        self.ttf          = None
        self.viseo        = 0.25

    def toString(self): 
        return self.mapper.toString()
    def register(self,mapper):
        self.mapper  = mapper
        self.mapper.register(self)
        pass
    def configure(self,**kwargs): 
        if self.mapping in ['t3d','oxal','base','ttf']:
            self.aperture  = kwargs.get('aperture')        # [m]
            self.cavlen    = kwargs.get('cavlen')          # [m] cavity length
            self.EzPeak    = kwargs.get('EzPeak')*self.dWf # [MV/m]
            self.freq      = kwargs.get('freq')            # [Hz]  RF frequenz
            self.gap       = kwargs.get('gap')             # [m] rf-gap
            self.phisoll   = kwargs.get('phisoll')         # [radians] soll phase
            self.sec       = kwargs.get('sec')
            self.SFdata    = kwargs.get('SFdata')          # SuperFish data
        
            self.lamb      = UTIL.PARAMS['clight']/self.freq
            self.omega     = twopi*self.freq
            # add attributes to kwargs and configure the mapper
            kwargs['lamb']     = self.lamb
            kwargs['omega']    = self.omega
            self.mapper.configure(**kwargs)
            pass
        elif self.mapping in ['simple','dyn']:
            raise(UserWarning(wrapRED(f'mapping not ready {self.mapping}')))
            sys.exit()
        else:
            raise(UserWarning(wrapRED(f'missing implementation {self.mapping}')))
            sys.exit()
    def adjust_energy(self, tkin):
        super().adjust_energy(tkin)
        self.mapper.adjust_energy(tkin)
        return self
    def map(self,i_track):
        return self.mapper.map(i_track)
    def make_slices(self, anz=1):
        return [self]
    def waccept(self):
        return self.mapper.waccept()
    def aper_check(self,new_tp,s,**kwargs):
        new_point=new_tp()            # calling bunch.Tpoint() object

        fifo_z   = kwargs['fifo_z']   # collect losses in tracker.Fifo objects
        sfifo_z  = kwargs['sfifo_z']
        fifo_xy  = kwargs['fifo_xy']
        sfifo_xy = kwargs['sfifo_xy']

        lost   = False
        mapper = self.mapper
        tkin   = mapper.particle.tkin

        res  = self.waccept()    # test energy acceptance
        (dummy,phi_2,phisoll,phi_1) = res['phaseacc']   # dummy kept for old track_node
        conv = UTIL.WConverter(tkin, mapper.freq)
        Dphi = conv.zToDphi(new_point[Ktp.z])
        phi  = mapper.phisoll+Dphi

        # longitudinal acceptance
        if not (phi_2 < phi and phi < phi_1):  # Wrangler's approximation (pp.178) is good up to -58 [deg]
            fifo_z.append(f'loss (z) {new_point[Ktp.z]:.3e} at {s:.4e} m')
            sfifo_z.append(s)
            lost = True
        elif abs(new_point[Ktp.zp]) > res['Dp2pmax']:
            fifo_z.append(f'loss (zp) {new_point[Ktp.zp]:.3e} at {s:.4e} m')
            sfifo_z.append(s)
            lost = True
        # transverse apertures
        elif mapper.aperture != None and not (abs(new_point[Ktp.x]) < mapper.aperture or abs(new_point[Ktp.y]) < mapper.aperture):
            fifo_xy.append(f'loss (x|y) ({new_point[Ktp.x]:.3e},{new_point[Ktp.y]:.3e}) at {s:.4e} m')
            sfifo_xy.append(s)
            lost = True
        return lost

class TestElementMethods(unittest.TestCase):
    def Matrix(self,v):
        matrix = NP.eye(MDIM,MDIM)
        matrix[0,0] = v[0]
        matrix[0,1] = v[1]
        matrix[1,0] = v[2]
        matrix[1,1] = v[3]
        matrix[2,2] = v[4]
        matrix[2,3] = v[5]
        matrix[3,2] = v[6]
        matrix[3,3] = v[7]
        matrix[4,4] = v[8]
        matrix[4,5] = v[9]
        matrix[5,4] = v[10]
        matrix[5,5] = v[11]
        return matrix
    def test_Node(self):
        """ testing the * operator for NODE objects"""
        print("\b----------------------------------------test_Node")
        a = self.Matrix((1,2,3,4,1,2,4,5,1,2,5,6))
        b = self.Matrix((1,0,0,1,1,0,0,1,1,0,0,1))
        b = 2*b
        c = self.Matrix((1,1,0,1,1,1,0,1,1,1,0,1))
        A = Node()
        A.matrix = a
        A.label = 'Node-A'
        A.length = 1.
        # print(A.toString())

        B = Node()
        B.matrix = b
        B.label = 'Node-B'
        B.length = 1.
        # print(B.toString())
        AB = A*B
        # print(AB.toString())
        BA = B*A
        # print(BA.toString())

        self.assertTrue(NP.array_equal(AB.matrix,a.dot(b)),msg='AB * ab')
        self.assertTrue(NP.array_equal(BA.matrix,b.dot(a)),msg='BA * ba')
        self.assertTrue(NP.array_equal(AB.matrix,BA.matrix),msg='A * B')

        C = Node()
        C.matrix = c
        C.label = 'Node-C'
        C.length = 1.
        # print(C.toString())
        AC = A*C
        # print(AC.toString())
        CA = C*A
        # print(CA.toString())

        self.assertFalse(NP.array_equal(AC.matrix,CA.matrix))
    def test_QF_Node_make_slices(self):
        print("\b----------------------------------------test_QF_Node_make_slices")
        gradient = 3.; length = 1.; p = UTIL.Proton(80.)
        QFnode = QF("QFoc", gradient, particle=p, length=length)
        # slices = QFnode.make_slices(anz=anz)
        anz=5
        self.assertEqual(len(QFnode.make_slices(anz=anz)),anz)
        for i in range(anz):
            self.assertEqual(QFnode.make_slices(anz=anz)[i].length,length/anz)
        anz=1
        self.assertEqual(len(QFnode.make_slices(anz=anz)),anz)
        for i in range(anz):
            self.assertEqual(QFnode.make_slices(anz=anz)[i].length,length/anz)
        anz=0
        self.assertEqual(len(QFnode.make_slices(anz=anz)),1)
        self.assertEqual(QFnode.make_slices(anz=anz)[0].length,length)
    def test_I_Node(self):
        print("\b----------------------------------------test_I_Node")
        i_matrix = NP.eye(MDIM,MDIM)
        Inode    = I("IN1")
        II       = Inode*Inode
        self.assertEqual(Inode.type,"I",'type')
        self.assertEqual(Inode.label,"IN1","label")
        self.assertEqual(Inode.length,0.,"length")
        self.assertTrue(NP.array_equal(Inode.matrix,i_matrix),"matrix")
        self.assertTrue(NP.array_equal(II.matrix,i_matrix.dot(i_matrix)),"I*I")
        self.assertEqual(II.label,"IN1*IN1")
        # NOTE: type( Child * Child ) = 'Node'
        self.assertTrue(II.type == "Node")
    def test_D_Node(self):
        print("\b----------------------------------------test_D_Node")
        l = 10.
        p = UTIL.Proton(50.)
        Dnode = D("Drift",length=l,aperture=5.)
        self.assertEqual(Dnode.length,10.)
        g = p.gamma
        d_matrix = NP.eye(MDIM,MDIM)
        d_matrix[XKOO,XPKOO] = d_matrix[YKOO,YPKOO] = d_matrix[SKOO,DSKOO] = l
        d_matrix[ZKOO,ZPKOO] = l/g**2
        self.assertTrue(NP.array_equal(Dnode.matrix,d_matrix),"matrix")
        tkin = 200.
        Dnode = Dnode.adjust_energy(tkin)
        p = p(tkin)
        g = p.gamma
        d_matrix[ZKOO,ZPKOO] = l/g**2
        self.assertTrue(NP.array_equal(Dnode.matrix,d_matrix),"adjust_energy")
        l = l/2.
        Dnode = Dnode.shorten(l)
        d_matrix[XKOO,XPKOO] = d_matrix[YKOO,YPKOO] = d_matrix[SKOO,DSKOO] = l
        d_matrix[ZKOO,ZPKOO] = l/g**2
        self.assertTrue(NP.array_equal(Dnode.matrix,d_matrix),"shorten")
        self.assertTrue(Dnode.label == "Drift")
        self.assertTrue(Dnode.type == "D")
        self.assertEqual(Dnode.particle.tkin,p.tkin,"tkin")
    def test_Particle_and_K(self):
        print("\b----------------------------------------test_Particle_and_K")
        p = UTIL.Proton(50.)
        gradient =3.
        K0 = K(gradient,p)
        self.assertAlmostEqual(p.tkin, 50., delta=1e-3)
        self.assertAlmostEqual(p.gamma, 1.0533, delta=1e-4)
        self.assertAlmostEqual(p.beta, 3.1405e-1, delta=1e-4)
        self.assertAlmostEqual(p.gamma_beta, 3.3079e-1, delta=1e-4)
        self.assertAlmostEqual(p.brho, 1.0353, delta=1e-4)
        self.assertAlmostEqual(K0, 2.8978, delta=1e-4)
    def test_QF_Node(self): 
        print("\b----------------------------------------test_QF_Node")
        def matrix(gradient,length,particle,thin):
            x=0; xp=1; y=2; yp=3; z=4; zp=5; T=6; dT=7; S=8; dS=9
            k02  = K(gradient,particle)
            k0w  = M.sqrt(k02)
            mx = NP.eye(MDIM,MDIM)
            if thin == False:
                dphi = k0w*length
                mx[x,x]  = M.cos(dphi);        mx[x,xp]  = M.sin(dphi)/k0w       
                mx[xp,x] = -k0w*M.sin(dphi);   mx[xp,xp] = mx[x,x]
                mx[y,y]  = M.cosh(dphi);       mx[y,yp]  = M.sinh(dphi)/k0w        
                mx[yp,y] = +k0w*M.sinh(dphi);  mx[yp,yp] = mx[y,y]
            else:
                l = length
                f = 1./(k02*l)
                mx[XKOO, XKOO]  = 1.-l/(2.*f); mx[XKOO, XPKOO] = l-l**2/(6*f); mx[XPKOO, XKOO] = -1./f; mx[XPKOO, XPKOO] = mx[XKOO, XKOO]
                f = -f
                mx[YKOO, YKOO]  = 1.-l/(2.*f); mx[YKOO, YPKOO] = l-l**2/(6*f); mx[YPKOO, YKOO] = -1./f; mx[YPKOO, YPKOO] = mx [YKOO, YKOO]
            mx[z,zp] = length/particle.gamma**2
            mx[S,dS] = length
            return mx
        gradient = 3.; length = 1.; p = UTIL.Proton(80.)
        QFnode = QF("QFoc", gradient, particle=p, length=length)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"matrix")
        self.assertEqual(QFnode.viseo,+0.5,"viseo")
        self.assertTrue(QFnode.matrix[XPKOO,XKOO] < 0.)
        self.assertTrue(QFnode.matrix[YPKOO,YKOO] > 0.)

        p = UTIL.Proton(100.)
        QFnode = QFnode.adjust_energy(100.)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"QF.adjust_energy")

        length = 0.5
        QFnode = QFnode.shorten(length)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"QF.shorten")

        gradient = 1.0; length = 0.15; p = UTIL.Proton(100.)
        QFnode = QF("QFfoc",gradient,particle=p,length=length)
        mx = matrix(gradient,length,p,QFnode.isThin())
        self.assertTrue(NP.array_equal(QFnode.matrix,mx),"matrix")
    def test_QD_Node(self): 
        print("\b----------------------------------------test_QD_Node")
        def matrix(gradient,length,particle,thin):
            x=0; xp=1; y=2; yp=3; z=4; zp=5; T=6; dT=7; S=8; dS=9
            k02  = K(gradient,particle)
            k0w  = M.sqrt(k02)
            mx = NP.eye(MDIM,MDIM)
            if thin == False:
                dphi = k0w*length
                mx[y,y]  = M.cos(dphi);        mx[y,yp]  = M.sin(dphi)/k0w        
                mx[yp,y] = -k0w*M.sin(dphi);   mx[yp,yp] = mx[y,y]
                mx[x,x]  = M.cosh(dphi);       mx[x,xp]  = M.sinh(dphi)/k0w        
                mx[xp,x] = +k0w*M.sinh(dphi);  mx[xp,xp] = mx[x,x]
            else:
                l = length
                f = -1./(k02*l)
                mx[XKOO, XKOO]  = 1.-l/(2.*f); mx[XKOO, XPKOO] = l-l**2/(6*f); mx[XPKOO, XKOO] = -1./f; mx[XPKOO, XPKOO] = mx[XKOO, XKOO]
                f = -f
                mx[YKOO, YKOO]  = 1.-l/(2.*f); mx[YKOO, YPKOO] = l-l**2/(6*f); mx[YPKOO, YKOO] = -1./f; mx[YPKOO, YPKOO] = mx [YKOO, YKOO]
            mx[z,zp] = length/particle.gamma**2
            mx[S,dS] = length
            return mx
        gradient = 3.; length = 1.; p = UTIL.Proton(80.)
        QDnode = QD("QDfoc", gradient, particle=p, length=length)
        mx = matrix(gradient,length,p,QDnode.isThin())
        self.assertTrue(NP.array_equal(QDnode.matrix,mx),"matrix")
        self.assertTrue(QDnode.matrix[XPKOO,XKOO] > 0.)
        self.assertTrue(QDnode.matrix[YPKOO,YKOO] < 0.)

        length = 0.5
        QDnode = QDnode.shorten(length)
        mx = matrix(gradient,length,p,QDnode.isThin())
        self.assertTrue(NP.array_equal(QDnode.matrix,mx),"QD.shorten")

        self.assertEqual(QDnode.label,'QDfoc',"__getitem__")
        self.assertEqual(QDnode.type,'QD',"__getitem__")
        self.assertEqual(QDnode.viseo,-0.5,"viseo")
        QDnode.viseo = 0.75
        self.assertEqual(QDnode.viseo,0.75,"__setitem__")

        gradient = 1.0; length = 0.15; p = UTIL.Proton(100.)
        QDnode = QD("QDfoc",gradient,particle=p,length=length)
        mx = matrix(gradient,length,p,QDnode.isThin())
        self.assertTrue(NP.array_equal(QDnode.matrix,mx),"matrix")
    def test_Wille_Nodes(self):
        print("\b----------------------------------------test_Wille_Nodes")
        kqf  = -1.20    # [1/m**2]
        kqd  = -kqf
        lqf  = 0.20       # [m]
        lqd  = 0.40       # [m]
        rhob = 3.8197     # [m]
        lb   = 1.50       #[m]
        phib = 11.25      # [deg]
        ld   = 0.55       # [m]
        p    = UTIL.Proton(100.)
        gradf = kqf*p.brho
        gradd = kqd*p.brho

        """ Wille's Fokussierende Quadrupole = MQF """
        MQF = QF("MQF",gradf,particle=p,length=lqf)
        mqf = NP.eye(MDIM,MDIM)
        mqf[XKOO,XKOO]  = 0.9761;     mqf[XKOO,XPKOO]  = 0.1984
        mqf[XPKOO,XKOO] = -0.2381;    mqf[XPKOO,XPKOO] = mqf[XKOO,XKOO]
        mqf[YKOO,YKOO]  = 1.0241;     mqf[YKOO,YPKOO]  = 0.2016
        mqf[YPKOO,YKOO] = 0.2419;     mqf[YPKOO,YPKOO] = mqf[YKOO,YKOO]
        mqf[ZKOO,ZPKOO] = lqf/p.gamma**2
        mqf[EKOO,EKOO]  = 1.;         mqf[EKOO,DEKOO]  = 0.     # tkin
        mqf[SKOO,SKOO]  = 1.;         mqf[SKOO,DSKOO]   = lqf     # s
        # print("\nmqf");print(mqf);print("QF.matrix");print(MQF.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mqf[i,j],MQF.matrix[i,j],msg="MQF",delta=1e-4)

        """ Wille's Defokussierende Quadrupole = MQD """
        MQD = QD("MQD",gradd,particle=p,length=lqd)
        mqd = NP.eye(MDIM,MDIM)
        mqd[XKOO,XKOO]  = 1.0975;     mqd[XKOO,XPKOO]  = 0.4129
        mqd[XPKOO,XKOO] = 0.4955;     mqd[XPKOO,XPKOO] = mqd[XKOO,XKOO]
        mqd[YKOO,YKOO]  = 0.9055;     mqd[YKOO,YPKOO]  = 0.3873
        mqd[YPKOO,YKOO] = -0.4648;    mqd[YPKOO,YPKOO] = mqd[YKOO,YKOO]
        mqd[ZKOO,ZPKOO] = lqd/p.gamma**2
        mqd[EKOO,EKOO]  = 1.;         mqd[EKOO,DEKOO]  = 0.     # tkin
        mqd[SKOO,SKOO]  = 1.;         mqd[SKOO,DSKOO]   = lqd     # s
        # print("\nmqd");print(mqd);print("MQD.matrix");print(MQD.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mqd[i,j],MQD.matrix[i,j],msg="MQD",delta=1e-4)

        """ Wille's Dipole = MB """
        MB = SD("MB",2*phib,rhob,particle=p)
        mb = NP.eye(MDIM,MDIM)
        mb[XKOO,XKOO]  = 0.9239;      mb[XKOO,XPKOO]   = 1.4617;        mb[XKOO,ZPKOO] = 0.2908
        mb[XPKOO,XKOO] = -0.1002;     mb[XPKOO,XPKOO]  = mb[XKOO,XKOO]; mb[XPKOO,ZPKOO] = 0.3827
        mb[YKOO,YKOO]  = 1.0;         mb[YKOO,YPKOO]   = 1.5
        mb[YPKOO,YKOO] = -0.0;        mb[YPKOO,YPKOO]  = mb[YKOO,YKOO]
        mb[ZKOO,XKOO]  = -0.3827;     mb[ZKOO,XPKOO]   = -0.2908;       mb[ZKOO,ZPKOO] = 1.1867
        mb[EKOO,EKOO]  = 1.;          mb[EKOO,DEKOO]   = 0.     # tkin
        mb[SKOO,SKOO]  = 1.;          mb[SKOO,DSKOO]    = lb     # s
        # print("\nmb");print(mb);print("MB.matrix");print(MB.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mb[i,j],MB.matrix[i,j],msg="MB",delta=1e-4)

        """ Wille's Kantenfokusierung = MEB """
        MEB = Wedge(phib,rhob, t3d_wedge=False)
        meb = NP.eye(MDIM,MDIM)
        meb[XPKOO,XKOO] = 0.0521; meb[YPKOO,YKOO] = -meb[XPKOO,XKOO]
        # print("\nmeb");print(meb);print("MEB.matrix");print(MEB.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(meb[i,j],MEB.matrix[i,j],msg="MEB",delta=1e-2)

        """ Wille's Driftstrecken = MD """
        MD = D("MD",particle=p,length=ld)
        md = NP.eye(MDIM,MDIM)
        md[XKOO,XPKOO] = ld; md[YKOO,YPKOO] = md[XKOO,XPKOO]
        md[ZKOO,ZPKOO] = ld/p.gamma**2
        md[SKOO,DSKOO] = ld
        # print("\nmd");print(md);print("MD.matrix");print(MD.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(md[i,j],MD.matrix[i,j],msg="MD",delta=1e-4)

        """ Wille's gesamte Zelle mit Sektordipol und Kantenfokussierung = MZ """
        MZ = MQF*MD*MEB*MB*MEB*MD*MQD*MD*MEB*MB*MEB*MD*MQF
        # MZ = MQF*MD*MEB*MB*MEB*MD*MQD   #*MD*MEB*MB*MEB*MD*MQF
        mz = NP.eye(MDIM,MDIM)
        mz[XKOO,XKOO]  = 0.0808;   mz[XKOO,XPKOO]  = 9.7855;        mz[XKOO,ZPKOO]  = 3.1424
        mz[XPKOO,XKOO] = -0.1015;  mz[XPKOO,XPKOO] = mz[XKOO,XKOO]; mz[XPKOO,ZPKOO] = 0.3471
        mz[YKOO,YKOO]  = -0.4114;  mz[YKOO,YPKOO]  = 1.1280
        mz[YPKOO,YKOO] = -0.7365;  mz[YPKOO,YPKOO] = mz[YKOO,YKOO]
        mz[ZKOO,XKOO]  = -0.34707; mz[ZKOO,XPKOO]  = -3.1424;       mz[ZKOO,ZPKOO]  = 4.1844
        mz[SKOO,DSKOO]  = 6.0
        # print("\nmz");print(mz);print("MZ.matrix");print(MZ.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mz[i,j],MZ.matrix[i,j],msg="MZ",delta=1e-4)
        """ Wille's gesamte Zelle mit RD-class """
        WEDGE = MEB
        MRD   = RD("RD",2*phib,rhob,WEDGE,particle=p)
        MRDZ  = MQF*MD*MRD*MD*MQD*MD*MRD*MD*MQF
        # print("\nmz");print(mz);print("MRDZ.matrix");print(MRDZ.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(mz[i,j],MRDZ.matrix[i,j],msg="MRDZ == MZ?",delta=1e-4)
    def test_SD_and_RD_Node_make_slices(self):
        print("\b----------------------------------------test_SD_and_RD_Node_make_slices")
        rhob   = 3.8197   # [m]
        phib   = 11.25    # [deg]
        p      = UTIL.Proton(100.)

        """ slice,recombine and adjust SD """
        sd     = SD("SD",2*phib,rhob,particle=p)
        slices = sd.make_slices(anz=3)
        # for item in slices: print(item.matrix); print()
        sdx = slices[0]         # combine the slices again
        for i in range(1,len(slices)):
            sdx = sdx * slices[i]
        # print("\nsd"); print(sd.matrix); print("\nsdx"); print(sdx.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(sd.matrix[i,j],sdx.matrix[i,j],msg='sd == sdx?',delta=1.e-4)

        sd_adjusted = sd.adjust_energy(tkin=p.tkin)
        for i in range(10):
            for j in range(10):
                self.assertEqual(sd.matrix[i,j],sd_adjusted.matrix[i,j],'sd == sd_adjusted?')

        """ slice,recombine and agjust RD """
        wedge =  Wedge(phib,rhob,t3d_wedge=True)
        rd = RD("RD",2*phib,rhob,wedge,particle=p)
        slices = rd.make_slices(anz=5)
        # for item in slices: print(item.matrix); print()
        rdx = slices[0]
        for i in range(1,len(slices)):
            rdx = rdx * slices[i]
        # print("\nrd"); print(rd.matrix); print("\nrdx"); print(rdx.matrix)
        for i in range(10):
            for j in range(10):
                self.assertAlmostEqual(rd.matrix[i,j],rdx.matrix[i,j],msg='rd == rdx?',delta=1.e-4)

        rd_adjusted = rd.adjust_energy(tkin=p.tkin)
        for i in range(10):
            for j in range(10):
                self.assertEqual(rd.matrix[i,j],rd_adjusted.matrix[i,j],'rd == rd_adjusted?')
    def test_SD_and_RD_Node_adjust_energy(self):
        print("\b----------------------------------------test_SD_and_RD_Node_adjust_energy")
        rhob   = 3.8197   # [m]
        phib   = 11.25    # [deg]
        p100      = UTIL.Proton(100.)
        p50       = UTIL.Proton(50.)

        """ adjust SD  100->50->100 """
        sd = SD("SD",2*phib,rhob,particle=p100)
        sd_adjusted = sd.adjust_energy(tkin=p50.tkin)
        sd_100 = sd_adjusted.adjust_energy(tkin=p100.tkin)
        for i in range(10):
            for j in range(10):
                self.assertEqual(sd.matrix[i,j],sd_100.matrix[i,j],'sd == sd_100?')

        """  adjust RD 100->50->100 """
        wedge =  Wedge(phib,rhob,t3d_wedge=True)
        rd = RD("RD",2*phib,rhob,wedge,particle=p100)
        rd_adjusted = rd.adjust_energy(tkin=p50.tkin)
        rd_100 = rd_adjusted.adjust_energy(tkin=p100.tkin)
        for i in range(10):
            for j in range(10):
                self.assertEqual(rd.matrix[i,j],rd_100.matrix[i,j],'rd == rd_100?')
    def test_GAP_Node(self):
        print("\b----------------------------------------test_GAP_Node")
        EzAvg   = 2.1             #[MV/m]
        phisoll = M.radians(-30.)
        gap     = 0.022           #[m]
        freq    = 816.e6          #[Hz]
        gap = GAP("GAP",EzAvg,phisoll,gap,freq,dWf=1) # tkin=500 default
        # print(gap.matrix)
        self.assertAlmostEqual(gap.matrix[EKOO,DEKOO],0.03766,delta=1.e-4)
        gap = gap.adjust_energy(6.) # tkin=6
        # print(gap.matrix)
        self.assertAlmostEqual(gap.matrix[EKOO,DEKOO],0.023817,delta=1.e-4)
    def test_RFG_Node_with_T3D_G_map(self):
        print("\b----------------------------------------test_RFG_Node_with_T3D_G_map")
        ID        = 'RFG mit T3D_G'
        EzPeak    = 2.1
        phisoll   = M.radians(-30.)
        gap_parameters = dict(
            EzPeak    = EzPeak,
            phisoll   = phisoll,         # [radians] requested soll phase
            gap       = 0.044,
            cavlen    = 0.064,
            freq      = 816.e6,            # [Hz]  requested RF frequenz
            particle  = Proton(50.),
            position  = (0,0,0),
            aperture  = 0.011,
            sec       = 'xx'
        )
        instance = RFG(ID)
        instance.register(T3D_G())
        instance.configure(**gap_parameters)
        
        self.assertEqual(instance.mapping,'t3d')
        self.assertEqual(instance.label,'RFG mit T3D_G')
        self.assertEqual(instance.mapper.label,'T3D_G')
        self.assertEqual(instance.mapper.EzPeak,EzPeak)
        self.assertEqual(instance.mapper.phisoll,M.radians(-30.))
        self.assertAlmostEqual(instance.mapper.deltaW,0.062206,delta=1.e-4)
        self.assertAlmostEqual(instance.mapper.matrix[EKOO,DEKOO],0.062206,delta=1.e-4)
        self.assertEqual(instance.length,0.)
        self.assertEqual(instance.mapper.matrix[SKOO,DSKOO],0.)
        instance.adjust_energy(250.)
        self.assertAlmostEqual(instance.mapper.deltaW,0.0751,delta=1.e-4)
        self.assertAlmostEqual(instance.matrix[EKOO,DEKOO],0.0751,delta=1.e-4)

if __name__ == '__main__':
    UTIL.FLAGS['verbose'] = 3
    UTIL.FLAGS['mapping'] = 't3d'
    unittest.main()



""" tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug  """
""" tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug tbd Zeug  """
# class SD(Node):
#     """ Trace3d horizontal sector magnet. n=0 pure dipole, alpha in [deg], rho in [m]."""
#     def __init__(self, label, alpha, rho, n=0, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None):
#         super().__init__()
#         self.label    = label
#         self.alpha    = alpha   # [deg]
#         self.rho      = rho     # [m]
#         self.n        = n
#         self._particle = copy(particle)
#         self.position = position
#         self.length   = rho*M.radians(self.alpha)
#         self.aperture = aperture
#         self.viseo    = 0.25
#         self.matrix   = self._mx()
#     def _mx(self):
#         m = NP.eye(MDIM,MDIM)
#         beta  = self.particle.beta
#         gamma = self.particle.gamma
#         alpha = M.radians(self.alpha)      # [rad]
#         rho   = self.rho                 # [m]
#         Ds     = abs(rho*alpha)          # [m]
#         h = alpha/Ds      # [1/m]
#         kx = M.sqrt((1-self.n)*h**2)       # [1/m]
#         ky = M.sqrt(self.n*h**2)
#         self.length = Ds
#         cx = M.cos(kx*Ds)
#         sx = M.sin(kx*Ds)
#         cy = M.cos(ky*Ds)
#         sy = M.sin(ky*Ds)

#         m[XKOO, XKOO]  = cx;      m[XKOO, XPKOO]    = sx/kx            # x,x'-plane
#         m[XPKOO, XKOO] = -kx*sx;  m[XPKOO, XPKOO]   = m[XKOO, XKOO]        
#         m[YKOO, YKOO]  = cy;      m[YKOO, YPKOO]    = sy/ky if self.n != 0. else self.length   # y,y'-plane
#         m[YPKOO, YKOO] = -ky*sy;  m[YPKOO, YPKOO]   = m[YKOO, YKOO]
#         # z,z'-plane
#         m[XKOO,ZKOO]  = 0.;       m[XKOO,ZPKOO]   = h*(1.-cx)/kx**2    # Rxz
#         m[XPKOO,ZKOO] = 0.;       m[XPKOO,ZPKOO]  = h*sx/kx
#         m[ZKOO,XKOO]  = -h*sx/kx; m[ZKOO,XPKOO]   = -h*(1.-cx)/kx**2   # Rzx
#         m[ZPKOO,XKOO] = 0.;       m[ZPKOO,XPKOO]  = 0.
#         m[ZKOO,ZKOO]  = 1.;       m[ZKOO,ZPKOO]   = -1./rho**2*(kx*Ds*beta**2-sx)/kx**3+Ds*(1.-1./(rho**2*kx**2))/gamma**2   #Rzz
#         m[ZPKOO,ZKOO] = 0.;       m[ZPKOO,ZPKOO]  = m[ZKOO,ZKOO]
#         m[SKOO, DSKOO]  = self.length  # length increase
#         return m
#     def adjust_energy(self, tkin):
#         adjusted = SD(self.label,self.alpha,self.rho,self.n,particle=self.particle(tkin),position=self.position,aperture=self.aperture)
#         return adjusted
#     def shorten(self, alpha):
#         shortSD = SD(self.label,alpha,self.rho,self.n, particle=self.particle,position=self.position, aperture=self.aperture)
#         return shortSD
#     def make_slices(self, anz=2):
#         shortSD = self.shorten(self.alpha/anz)
#         slices = []
#         for i in range(anz):
#             slices.append(shortSD)
#         return slices
# class RD(SD):
#     """ Trace3d horizontal rechteck magnet. n=0 pure dipole, alpha in [deg], rho in [m]."""
#     def __init__(self, label, alpha, rho, wedge, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None):
#         super().__init__(label, alpha, rho, n=0, particle=particle, position=position, aperture=aperture)
#         self.wedge  = wedge
#         self.matrix = NP.dot(self.wedge.matrix,NP.dot(self.matrix,self.wedge.matrix))
#     def adjust_energy(self, tkin):
#         adjusted = RD(self.label,self.alpha,self.rho,self.wedge,particle=self.particle(tkin),position=self.position,aperture=self.aperture)
#         return adjusted
#     def shorten(self, alpha):
#         shortSD = RD(self.label,alpha,self.rho,self.wedge, particle=self.particle,position=self.position, aperture=self.aperture)
#         return shortSD
#     def make_slices(self, anz=2):
#         if anz < 2: anz = 2
#         slicewinkel   = self.alpha/anz
#         sdi   = SD(self.label,slicewinkel,self.rho,particle=self.particle,position=self.position,aperture=self.aperture)
#         sdi.matrix = NP.dot(self.wedge.matrix,sdi.matrix)
#         slices  = [sdi]          # wedge @ entrance
#         for i in range(1,anz-1):
#             shortRD = SD(self.label,slicewinkel,self.rho,particle=self.particle,position=self.position,aperture=self.aperture)
#             slices.append(shortRD)
#         sdf   = SD(self.label,slicewinkel,self.rho,particle=self.particle,position=self.position,aperture=self.aperture)
#         sdf.matrix = NP.dot(sdf.matrix,self.wedge.matrix)
#         slices.append(sdf)
#         return slices
# class Wedge(Node):
#     """  Trace3d kantenfokussierung .a.k.a wedge: Ryy simplified if t3d_wedge=False """
#     def __init__(self, kwinkel, rho, t3d_wedge=True):
#         super().__init__()
#         self.kwinkel = kwinkel     # [deg ] kantenwinkel
#         self.rho = rho
#         self.t3d_wedge = t3d_wedge
#         self.label = "W"
#         self.length = 0.
#         beta = M.radians(self.kwinkel)    
#         g = 0.050    # fixed gap of 50 mm assumed
#         K1=0.45
#         K2=2.8
#         psi  = K1*g/rho*((1+M.sin(beta)**2)/M.cos(beta))*(1-K1*K2*(g/rho)*M.tan(beta)) if self.t3d_wedge else 0.
#         mxpx = M.tan(beta)/rho
#         mypy = M.tan(beta-psi)/rho
#         mx   = NP.eye(MDIM,MDIM)
#         mx[XPKOO, XKOO] = mxpx
#         mx[YPKOO, YKOO] = -mypy
#         self.matrix = mx
# class GAP(Node):
    # """ Simple zero length RF-gap nach Dr.Tiede & T.Wrangler
    # ... nicht sehr nuetzlich: produziert keine long. Dynamik wie Trace3D RFG!  """
    # def __init__(self, label, EzAvg, phisoll, gap, freq, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=UTIL.FLAGS['dWf']):
    #     """ EzAvg [MV/m], phisoll [rad], gap [m], freq [Hz] """
    #     super().__init__()
    #     self.accelerating = True
    #     self.label    = label
    #     self.EzAvg    = EzAvg*dWf
    #     self.phisoll  = phisoll
    #     self.gap      = gap
    #     self.freq     = freq
    #     self._particle = copy(particle)
    #     self.position = position
    #     self.aperture = aperture
    #     self.dWf      = dWf
    #     self.viseo    = 0.25
    #     self.length   = 0.
    #     E0L           = self.EzAvg*self.gap            # Spaltspannung
    #     self.deltaW,self.matrix = self._mx(E0L)
    # def _mx(self,E0L):
    #     """ cavity nach Dr.Tiede pp.33 """
    #     lamb   = UTIL.PARAMS['clight']/self.freq            # [m] wellenlaenge
    #     beta   = self.particle.beta                    # beta Einstein
    #     ttf    = self.ttf(beta,lamb,self.gap)          # time-transition factor
    #     bg     = self.particle.gamma_beta
    #     m0c2   = self.particle.e0
    #     deltaW = E0L*ttf*M.cos(self.phisoll)                   # delta-W T.Wrangler pp.221
    #     m      = NP.eye(MDIM,MDIM)
    #     cyp = cxp = -UTIL.pi*E0L*ttf*M.sin(self.phisoll)/(m0c2*lamb*bg**3)
    #     m[XPKOO, XKOO] = cxp
    #     m[YPKOO, YKOO] = cyp
    #     m[EKOO, DEKOO] = deltaW      # energy increase
    #     m[SKOO, DSKOO]  = self.length # length increase
    #     return deltaW,m
    # def ttf(self, beta, lamb, gap):
    #     """ ttf-factor nach Panofsky (Lapostolle CERN-97-09 pp.65, T.Wangler pp.39) """
    #     x = UTIL.pi*gap/(beta*lamb)
    #     ttf = NP.sinc(x/UTIL.pi)
    #     return ttf
    # def adjust_energy(self, tkin):
    #     adjusted = GAP(self.label,self.EzAvg,self.phisoll,self.gap,self.freq,particle=self.particle(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf)
    #     return adjusted
