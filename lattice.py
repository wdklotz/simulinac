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
# import traceback
from math import sqrt,fabs,acos,asin,pi,degrees
from numpy import linalg as LA
import numpy as NP
from copy import copy
import warnings

from setutil import wille,PARAMS,FLAGS,SUMMARY,objprnt,printv,DEBUG,mxprnt,KeepValues
from elements import XKOO,XPKOO,YKOO,YPKOO,ZKOO,ZPKOO,EKOO,DEKOO,SKOO,LKOO
import elements as ELM
from sigma import Sigma
import TTFG as TTF

## DEBUG MODULE
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF

## Lattice
class Lattice(object):
    """
    The Lattice object is a list of tuples: 
        (ELM.<element>, from_position, to_position)
    """
    def __init__(self):
        self.seq    = []
        self.length = 0.
        self.betax0 = 0.
        self.alfax0 = 0.
        self.gammx0 = 0.
        self.betay0 = 0.
        self.alfay0 = 0.
        self.gammy0 = 0.
        self.accel  = 0.

    def add_element(self,elment):
        """
        Add element to lattice
        """
        if len(self.seq) == 0:
            s0 = 0.
        else:
            s0 = self.seq[-1][-1]
        l = elment.length
        self.length += l
        elm_with_position = (elment,s0,s0+l)
        self.seq.append(elm_with_position)

    def string(self):
        """
        Log lattice layout to string (could be even better?)
        """
        mcell = ELM.I(label='')   ##  chain matrices
        for ipos in self.seq:
            element,s0,s1 = ipos
            printv(3,'{:10s}({:d})\tlength={:.3f}\tfrom-to: {:.3f} - {:.3f}'.
                format(element.label,id(element),element.length,s0,s1))
            mcell = element * mcell   ## Achtung: Reihenfolge im Produkt ist wichtig! Umgekehrt == Blödsinn
            mcell.set_section('<= full lattice map')
        return mcell.string()

    def stats(self,soll_track):
        """
        Gather lattice statistics
        """
        cav_counter = 0
        qf_counter  = 0
        qd_counter  = 0
        ttfm = +1.e+50
        ttfx = +1.e-50
        tk_i = soll_track.first()[6]
        tk_f = soll_track.last()[6]
        for item in self.seq:
            element,s0,s1 = item
            if isinstance(element,ELM.QF) and (not isinstance(element,ELM.QD)):
                qf_counter += 1
            if isinstance(element,ELM.QD):
                qd_counter += 1
            if isinstance(element,ELM.RFG) \
                or isinstance(element,ELM.RFC) \
                or isinstance(element,TTF.TTFG):
                cav_counter += 1
                ttfm = min(element.tr,ttfm)
                ttfx = max(element.tr,ttfx)
        SUMMARY['nbof F-quads*']        = qf_counter
        SUMMARY['nbof D-quads*']        = qd_counter
        SUMMARY['nbof cavities*']       = cav_counter
        SUMMARY['(ttf)min,(ttf)max*']   = (ttfm,ttfx)
        SUMMARY['(energy)i,(energy)f [MeV]']  = (tk_i,tk_f)

    def cell(self,closed=True):
        """
        Construct the full lattice-cell matrix and extract standard quantities:
            full cell: mcell
            stability?
            betatron tunes: mux, muy 
            det(M)
            check symplecticity
            twiss prameters beta, alpha, gamma for periodic lattices
            
        """
        mcell = ELM.I(label=' <==')   ##  chain matrices for full cell
        for count, ipos in enumerate(self.seq):
            element,s0,s1 = ipos
            mcell = element * mcell   ## Achtung: Reihenfolge im Produkt ist wichtig! Umgekehrt == Blödsinn

        ## Stabilität ?
        unstable = False
        stab = fabs(mcell.tracex())
        # if verbose:
        printv(1,'stability X? ',stab)
        if stab >= 2.0:
            # if verbose:
            printv(1,'unstable Lattice in x-plane\n')
            unstable = True
        else:
            cos_mux = 0.5 * stab
            mux = degrees(acos(cos_mux))

        stab = fabs(mcell.tracey())
        # if verbose:
        printv(1,'stability Y? ',stab)
        if stab >= 2.0:
            # if verbose:
            printv(1,'unstable Lattice in y-plane\n')
            unstable = True
        else:
            cos_muy = 0.5 * stab
            muy = degrees(acos(cos_muy))
        if not unstable:
            # if verbose:
            printv(1,'\nphase_advance: X[deg]={:3f} Y[deg]={:.3f}\n'.format(mux,muy))
        ## full accelerator
        self.accel = mcell    # the full cell becomes instance variable
        # if verbose:
        printv(0,'Full Accelerator Matrix (f)<==(i)')
        printv(0,self.accel.string())
        det = LA.det(self.accel.matrix)
        # if verbose:
        printv(2,'det|full-cell|={:.5f}\n'.format(det))
        ## Determinate M-I == 0 ?
        beta_matrix = mcell.beta_matrix()
        for i in range(5):
            beta_matrix[i,i] = beta_matrix[i,i]-1.0
        det = LA.det(beta_matrix)
        # if verbose:
        printv(2,'det|Mbeta - I|={:.5f}\n'.format(det))
        ## symplectic?
        s = self.symplecticity()
        # if verbose:
        printv(2,'symplectic (+1,-1,+1,-1,+1,-1)?')
        printv(2,'[{:4>+.2f}, {:4>+.2f}, {:4>+.2f}, {:4>+.2f}, {:4>+.2f}, {:4>+.2f}]\n'.
            format(s[0],s[1],s[2],s[3],s[4],s[5]))
        ## Vorgabe emittance @ entrance
        emix = PARAMS['emitx_i']  
        emiy = PARAMS['emity_i']

        ## Startwerte für twiss-functions aus cell matrix (not beta_matrix!)
        if closed:
            if not unstable:
                cell_matrix = self.accel.matrix
                # m11  = cell_matrix[0,0];         m12  = cell_matrix[0,1]
                # m21  = cell_matrix[1,0];         m22  = cell_matrix[1,1]
                # n11  = cell_matrix[2,2];         n12  = cell_matrix[2,3]
                # n21  = cell_matrix[3,2];         n22  = cell_matrix[3,3]
                m11  = cell_matrix[XKOO,XKOO];   m12  = cell_matrix[XKOO,XPKOO]
                m21  = cell_matrix[XPKOO,XKOO];  m22  = cell_matrix[XPKOO,XPKOO]
                n11  = cell_matrix[YKOO,YKOO];   n12  = cell_matrix[YKOO,YPKOO]
                n21  = cell_matrix[YPKOO,YKOO];  n22  = cell_matrix[YPKOO,YPKOO]
                ## betax,alphax,gammax from transfer matrix 
                #  [m11,m12] = [cos(mu) + alpha * sin(mu),       beta * sin(mu)     ]
                #  [m21,m22] = [     -gamma * sin(mu)    , cos(mu) - alpha * sin(mu)]
                #  und: beta * gamma - alpha^2 = 1
                sin2mu = -(((m11-m22)**2)/4.+m12*m21)
                # print('sin(mux)^2={:4.4f}'.format(sin2mu))
                sinmu = sqrt(sin2mu)
                if m12 < 0:     # get rid of +- ambiguity from sqrt(sin(mu)^2)
                    sinmu = -sinmu                    
                bax = m12/sinmu 
                gmx = -m21/sinmu 
                alx2 = bax*gmx-1.
                # print('alfax^2',alx2)
                alx = 0. if fabs(alx2) < 1.e-9 else sqrt(alx2)
                print('betax {:4.4f} alfax {:4.4f} gammax {:4.4f}'.format(bax,alx,gmx))
                ## betay,alphay,gammay from transfer matrix                 
                sin2mu = -((n11-n22)**2/4.+n12*n21)
                # print('sin(muy)^2={:4.4f}'.format(sin2mu))
                sinmu = sqrt(sin2mu)
                if n12 < 0:     # get rid of +- ambiguity from sqrt(sin(mu)^2)
                    sinmu = -sinmu                    
                bay = n12/sinmu
                gmy = -n21/sinmu
                aly2 = bay*gmy-1.
                # print('alfay^2',aly2)
                aly = 0. if fabs(aly2) < 1.e-9 else sqrt(aly2)
                print('betay {:4.4f} alfay {:4.4f} gammay {:4.4f}'.format(bay,aly,gmy))                
                ## Probe: twiss-functions durch ganze Zelle mit beta-matrix (nur sinnvoll fuer period. Struktur!)
                v_beta_a = NP.array([bax,alx,gmx,bay,aly,gmy])
                m_cell_beta = self.accel.beta_matrix()
                v_beta_e = m_cell_beta.dot(v_beta_a)
                # if verbose:
                printv(1,'Probe: {TW(f)} == {BetaMatrix}x{TW(i)}?')
                diffa_e = v_beta_a - v_beta_e
                for i in range(6):
                    if fabs(diffa_e[i]) < 1.e-9: diffa_e[i] = 0.
                printv(1,'TW(i)-TW(f) (should be [0,...,0]):\n',diffa_e)
                ## keep variables for later use
                PARAMS['sigx_i'] = sqrt(bax*emix)
                PARAMS['sigy_i'] = sqrt(bay*emiy)
                SUMMARY['(sigx)i [mm]'] = 1000.*PARAMS['sigx_i']
                SUMMARY['(sigy)i [mm]'] = 1000.*PARAMS['sigy_i']
                xip = sqrt(emix*gmx)   # 1 sigma x' particle divergence
                yip = sqrt(emiy*gmy)
                SUMMARY["(sigx')i* [mrad]"] = 1000.*xip
                SUMMARY["(sigy')i* [mrad]"] = 1000.*yip
            else:
                print('unstable lattice - STOP')
                sys.exit(1)
        else:
            # Startwerte fuer transfer line (keine periodischen Randbedingungen!)
            # alfa, beta und emittance definieren den beam @ entrance
            # NOTE: transfer lattices need not to be stable!
            bax = PARAMS['betax_i']  # twiss beta @ entrance
            bay = PARAMS['betay_i']
            alx = PARAMS["alfax_i"]  # twiss alpha @ entrance
            aly = PARAMS["alfay_i"]
            gmx = (1.+alx*alx)/bax  # twiss gamma @ entrance
            gmy = (1.+aly*aly)/bay
            xip = sqrt(emix*gmx)   # 1 sigma x' particle divergence @ entrance
            yip = sqrt(emiy*gmy)
            SUMMARY["(sigx')i* [mrad]"] = 1000.*xip
            SUMMARY["(sigy')i* [mrad]"] = 1000.*yip
        ## keep twiss values as lattice instance varibles
        self.betax0 = bax
        self.alfax0 = alx
        self.gammx0 = gmx
        self.betay0 = bay
        self.alfay0 = aly
        self.gammy0 = gmy
        printv(0,'using @ entrance: [beta,  alfa,  gamma]-X    [beta,   alfa,   gamma]-Y')
        printv(0,'                  [{:.3f}, {:.3f}, {:.3f}]-X    [{:.3f},  {:.3f},  {:.3f}]-Y'.format(bax,alx,gmx,bay,aly,gmy))
        return (self.accel,self.betax0,self.betay0)

    def report(self):
        """
        Report lattice layout (may not work!)
        """
        raise RuntimeWarning('Lattice.report() not ready')
        reprt = ''
        header = ''
        row = ''
        for count, ipos in enumerate(reversed(self.seq)):
            elm,si,sf = ipos
            name = elm.label
            len = elm.length
            rest = (count+1)%19
            if rest != 0:
                header += '{:6s}'.format(name)
                row += '{:.3f} '.format(len)
            else:
                header += '{:6s}\n'.format(name)
                row += '{:.3f} \n'.format(len)
                reprt += header+row
                header = row = ''
        reprt += header+' \n'+row+' \n'
        return reprt

    def reverse(self):
        raise RuntimeWarning('Lattice.reverse() not implemented and not used! (probably bogus!)')
        res = Lattice()
        seq = copy(self.seq)
        seq.reverse()
        for ipos in seq:
            elm,s,s = ipos
            res.add_element(elm)
        return res

    # def append(self,piece):
    def concat(self,piece):
        """
        Concatenate two Lattice pieces
        """
        seq = copy(piece.seq)
        for ipos in seq:
            elm,s0,s1 = ipos
            self.add_element(elm)

    def twiss_functions(self,steps=10):
        """
        Track twiss functions with beta-matrix through lattice and scale to sigmas
        """
        beta_fun = []
        bx = self.betax0
        ax = self.alfax0
        gx = self.gammx0
        by = self.betay0
        ay = self.alfay0
        gy = self.gammy0
        v_beta0 = NP.array([bx,ax,gx,by,ay,gy])
        v_beta = v_beta0
        s = 0.0
        for ipos in self.seq:
            element,s0,s1 = ipos
            # particle = element.particle                                      # DEBUG
            # objprnt(particle,text='twiss_functions: '+element.label)         # DEBUG
            slices = element.make_slices(anz=steps)
            for i_element in slices:
                m_beta = i_element.beta_matrix()
                v_beta = m_beta.dot(v_beta)
                s += i_element.length
                betax  = v_beta[0]
                betay  = v_beta[3]
                viseo  = i_element.viseo
                beta_fun.append((s,betax,betay,viseo))
        return beta_fun

    def sigma_functions(self,steps=10):
        """
        Track the sigma-matrix through the lattice
        """
        sigma_fun = []
        # sigma initial
        sigma_i = Sigma(emitx=PARAMS['emitx_i'], betax=self.betax0,    alphax=self.alfax0,
                        emity=PARAMS['emity_i'], betay=self.betay0,    alphay=self.alfay0,
                        emitz=PARAMS['emitz_i'], betaz=PARAMS['betaz_i'],alphaz=0.)
        saper = 1.e6  # aperture control
        s     = 0.0
        for ipos in self.seq:
            element,s0,s1 = ipos
            # objprnt(element.particle ,text='sigma_functions: '+element.label)         # DEBUG
            slices = element.make_slices(anz=steps)
            for i_element in slices:
                # DEBUG_MODULE('{} {} {}'.format(i_element.__class__.__name__,'s0,s1',(s0,s1)))
                sigma_f = sigma_i.RSRt(i_element)        # map: sigma_f = R*sigma_i*RT
                if isinstance(i_element,ELM.RFG) and FLAGS['egf']:
                    rf_gap    = i_element
                    delta_phi = PARAMS['Dphi0']
                    sigma_f   = sigma_f.apply_eg_corr(rf_gap,sigma_i,delta_phi)
                sigf = sigma_f.matrix
                try:
                    xsquare_av = sqrt(sigf[0,0])   # sigmax = <x*x>**1/2 [m]
                    ysquare_av = sqrt(sigf[2,2])   # sigmay = <y*y>**1/2 [m]
                except ValueError:
                    # traceback.format_exc(limit=2)
                    # print('WARNING: results may not have physical meanings!')
                    warnings.showwarning(
                            'WARNING: sqrt of negative number!\nResult may have no physical meaning!',
                            UserWarning,
                            'lattice.py',
                            'sigma_functions()',
                            # line="xsquare_av = sqrt(sigf[0,0])   # sigmax = <x*x>**1/2 [m]"
                            )
                r = sqrt(xsquare_av**2+ysquare_av**2)
                s += i_element.length
                viseo = i_element.viseo
                sigma_fun.append((s,xsquare_av,ysquare_av,viseo))
                KeepValues.update({'z':s,'sigma_x':xsquare_av,'sigma_y':ysquare_av,'Tkin':i_element.particle.tkin})   # keep current values
                sigma_i = sigma_f.clone()
                if isinstance(i_element,ELM.MRK):                        # marker actions
                    i_element.do_actions()
            if 3.*r > PARAMS['aperture']:    # aperture control
                saper = min(s0,saper)
        if saper<1.e6:        # make use of warnings (experimental!)
            warnings.showwarning(
                    '3*sigma out of APERTURE at about s ={:5.1f}[m]\nParticle lost!'.format(saper),
                    UserWarning,
                    'lattice.py',
                    'sigma_functions()',
                    # line="if 3.*r > PARAMS['aperture']:    # aperture control"
                    )
        return sigma_fun

    def dispersion(self,steps=10,closed=True): 
        """
        Track the dispersion function
        """
        traj = []
        v_0 = NP.array([0.,0.,0.,0.,0.,1.,0.,0.,0.,0.])
        v_0.shape = (ELM.MDIM,1)   # column vector with MDIM rows, 1 column
        if closed == True:
            m_cell = self.accel
            m11 = m_cell.matrix[0,0]
            m15 = m_cell.matrix[0,5]
            d0  =  m15/(1.-m11)     # from H.Wiedemann (6.79) pp.206
            v_0[0,0] = d0
        s = 0.0
        for ipos in self.seq:
            element,s0,s1 = ipos
            slices = element.make_slices(anz=steps)
            for i_element in slices:
                m_beta = i_element.matrix
                v_0 = m_beta.dot(v_0)
                s += i_element.length
                d  = v_0[0,0]
                dp = v_0[1,0]
                viseo = i_element.viseo
                traj.append((s,d,dp))
        return traj

    def cs_traj(self,steps=10):
        """
        Track Cos & Sin trajectories
        """
        def SollTest_ON(arg):  # set all 0. to simulate Sollteilchen
            return 0.
        def SollTest_OFF(arg):
            return arg        
        soll_test   = SollTest_OFF

        gamma       = PARAMS['sollteilchen'].gamma
        beta        = PARAMS['sollteilchen'].beta
        tkin        = PARAMS['sollteilchen'].tkin
        lamb        = PARAMS['wellenlänge']
        x1          = soll_test(sqrt(PARAMS['emitx_i']*self.betax0)) # x-plane: principal-1 (cos like)
        x2p         = soll_test(sqrt(PARAMS['emitx_i']*self.gammx0)) # x-plane: principal-1 (sin like)
        y1          = soll_test(sqrt(PARAMS['emity_i']*self.betay0))
        y2p         = soll_test(sqrt(PARAMS['emity_i']*self.gammy0))
        sigmaz_i    = soll_test(PARAMS['sigmaz_i'])                  # z[m] Vorgabe
        dp2p_i      = soll_test(PARAMS['dp2p_i']*1.e-2)              # dp/p[%] Vorgabe --> []
        # MDIM tracking used here
        c_like = []
        s_like = []
        c_0 = NP.zeros(ELM.MDIM)
        s_0 = NP.zeros(ELM.MDIM)
        c_0[XKOO]  = x1; c_0[YKOO]  = y1;  c_0[ZKOO]  = sigmaz_i; c_0[EKOO] = tkin; c_0[DEKOO] = 1.; c_0[LKOO] = 1.  # cos-like traj.
        s_0[XPKOO] =x2p; s_0[YPKOO] = y2p; s_0[ZPKOO] = dp2p_i  ; s_0[EKOO] = tkin; s_0[DEKOO] = 1.; s_0[LKOO] = 1.  # sin-like traj.
        for ipos in self.seq:
            element,s0,s1 = ipos
            particle = element.particle
            gamma = particle.gamma
            # objprnt(particle,text='cs_traj: '+element.label)         # DEBUG
            slices = element.make_slices(anz=steps)
            for i_element in slices:
                ## cos_like
                # DEBUG_MODULE('cs_traj calls {}.map() for C'.format(i_element))
                c_0 = i_element.map(c_0)
                cx  = c_0[XKOO]
                cxp = c_0[XPKOO]
                cy  = c_0[YKOO]
                cyp = c_0[YPKOO]
                cz  = -c_0[ZKOO]*360./(beta*lamb)            # conversion sigmaz_i --> dPhi [deg]
                cdw = c_0[ZPKOO]*(gamma+1.)/gamma*100.       # conversion dp/p --> dW/W [%]
                ## sin_like
                # DEBUG_MODULE('cs_traj calls {}.map() for S'.format(i_element))
                s_0 = i_element.map(s_0)
                sx  = s_0[XKOO]
                sxp = s_0[XPKOO]
                sy  = s_0[YKOO]
                syp = s_0[YPKOO]
                sz  = -s_0[ZKOO]*360./(beta*lamb)
                sdw = s_0[ZPKOO]*(gamma+1.)/gamma*100.
                c_like.append((cx,cxp,cy,cyp,cz,cdw))
                s_like.append((sx,sxp,sy,syp,sz,sdw))
        return (c_like,s_like)

    def symplecticity(self):
        """
        Test symplecticity
        """
        s = NP.array([
                    [ 0.,1., 0.,0., 0.,0.,0.,0.,0.,0.],    #x
                    [-1.,0., 0.,0., 0.,0.,0.,0.,0.,0.],    #x'
                    [ 0.,0., 0.,1., 0.,0.,0.,0.,0.,0.],    #y
                    [ 0.,0.,-1.,0., 0.,0.,0.,0.,0.,0.],    #y'
                    [ 0.,0., 0.,0., 0.,1.,0.,0.,0.,0.],    #z
                    [ 0.,0., 0.,0.,-1.,0.,0.,0.,0.,0.],    #z'
                    [ 0.,0., 0.,0., 0.,0.,1.,0.,0.,0.],    #delta-E
                    [ 0.,0., 0.,0., 0.,0.,0.,1.,0.,0.],    #1
                    [ 0.,0., 0.,0., 0.,0.,0.,0.,1.,0.],    #delta-l
                    [ 0.,0., 0.,0., 0.,0.,0.,0.,0.,1.]     #1
                    ])
        s = NP.dot(self.accel.matrix.T,s)
        s = NP.dot(s,self.accel.matrix)
        # dets = LA.det(s)
        # if fabs(dets-1.) > 1.e-12:
            # for i in range(ELM.Matrix._dim):
                # print('[{:4>+.6f}, {:4>+.6f}, {:4>+.6f}, {:4>+.6f}, {:4>+.6f}, {:4>+.6f}]\n'.
                    # format(s[i,0],s[i,1],s[i,2],s[i,3],s[i,4],s[i,5]),end='')
        res = [s[0,1],s[1,0],s[2,3],s[3,2],s[4,5],s[5,4]]
        return(res)
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
## utilities
def make_wille():  
    """
    Wille's test lattice
    """
    print("K.Wille's Beispiel auf pp. 112-113")
    kqf = wille()['k_quad_f']
    lqf = wille()['length_quad_f']
    kqd = wille()['k_quad_d']
    lqd = wille()['length_quad_d']
    rhob = wille()['bending_radius']
    lb = wille()['dipole_length']
    ld = wille()['drift_length']
    ## elements
    mqf = ELM.QF(kqf,lqf,'QF')
    mqd = ELM.QD(kqd,lqd,'QD')
    mb  = ELM.SD(rhob,lb,'B')
    mb1 = ELM.SD(rhob,lb*0.5,'B1')  ## 1/2 sector dip.
    mw  = ELM.WD(mb)
    mw1 = ELM.WD(mb1)
    mbr = ELM.RD(rhob,lb)
    md  = ELM.D(ld)
    ## lattice
    lattice = Lattice()
    lattice.add_element(mqf)
    lattice.add_element(md)
    # lattice.add_element(mw)
    # lattice.add_element(mb)
    # lattice.add_element(mw)
    lattice.add_element(mbr)
    # lattice.add_element(mw1)
    # lattice.add_element(mb1)
    # lattice.add_element(mw1)
    lattice.add_element(md)
    lattice.add_element(mqd)
    lattice.add_element(md)
    # lattice.add_element(mw)
    # lattice.add_element(mb)
    # lattice.add_element(mw)
    lattice.add_element(mbr)
    lattice.add_element(md)
    lattice.add_element(mqf)
    # lattice.string()
    top = Lattice()
    top.concat(lattice)
    # top.concat(top)
    # top.concat(top)
    # top.concat(top)
    # top.concat(top)
    # top.concat(top)
    # top.concat(top)
    return top

def test0():
    print('-------------------------------------Test0--')
    lat = make_wille()
    print(lat.string())

def test1():
    from matplotlib.pyplot import plot,show,legend
    print('-------------------------------------Test2--')
    lattice = make_wille()
    # cell boundaries
    mcell,betax,betay = lattice.cell(closed=True)
    lattice.symplecticity()
    # twiss functions
    beta_fun = lattice.twiss_functions(steps=100)
    # cl,sl = lattice.cs_traj(steps=100)sK
    disp = lattice.dispersion(steps=100,closed=True)
    # plots
    vsbase = -1.
    s  = [x[0] for x in beta_fun]    # s
    xs = [x[1] for x in beta_fun]    # betax
    ys = [x[2] for x in beta_fun]    # betay
    ds = [x[1] for x in disp]        # dispersion
    vs = [x[3]+vsbase for x in beta_fun]  # viseo
    zero = [vsbase for x in beta_fun]     # viseo base line

    plot(s,xs,label='bx/bx0')
    plot(s,ys,label='by/by0')
    plot(s,ds,label='dp/p')
    plot(s,vs,label='element',color='black')
    plot(s,zero,color='black')
    legend(loc='upper left')
    show()
## main ----------
if __name__ == '__main__':
    test0()
    test1()
