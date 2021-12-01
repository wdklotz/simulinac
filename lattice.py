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
import numpy as NP
import warnings
import pprint, inspect
from math import sqrt,fabs,acos,degrees
from numpy import linalg as LA
from copy import copy

def PRINT_PRETTY(obj=None):
    file = inspect.stack()[0].filename
    print('DEBUG_ON ==============>  '+file)
    if obj != None: pprint.PrettyPrinter(width=200,compact=True).pprint(obj)
def PASS(obj=None):
    pass
DEB = dict(OFF=PASS,ON=PRINT_PRETTY)
DEBUG_ON  = DEB.get('ON')
DEBUG_OFF = DEB.get('OFF')

import elements as ELM
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO
from setutil import wille,PARAMS,FLAGS,SUMMARY,printv,sigmas, objprnt, Ktw, Ktp
from setutil import Twiss, Functions
from sigma import Sigma
# Lattice
class Lattice(object):
    """
    The Lattice object is a list of elements: ELM.<element> in self.seq  ??
    """
    class LRiterator(object):
        def __init__(self,lattice):
            self.lattice = lattice
            self.next = None
            if len(self.lattice.seq) > 0:
                self.next = self.lattice.seq[0]
        def __iter__(self):
            return self
        def __next__(self):
            if self.next != None:
                this = self.next
                self.next = self.next.next
                return this
            else:
                raise StopIteration
    class RLiterator(object):
        def __init__(self,lattice):
            self.lattice = lattice
            self.next = None
            if len(self.lattice.seq) > 0:
                self.next = self.lattice.seq[-1]
        def __iter__(self):
            return self
        def __next__(self):
            if self.next != None:
                this = self.next
                self.next = self.next.prev
                return this
            else:
                raise StopIteration    
    def __init__(self):
        self.seq    = []       # list of _Node objects z.B. [D,QD,GAP,QF....]
        self.iteration = "LR"  # default: iterating lattice left-right
        self.length = 0.
        self.accel  = 0.
    def __iter__(self):
        """ iterator using the linked list of element """
        if self.iteration == "RL":
            iterator = self.RLiterator(self)
        elif self.iteration == "LR":
            iterator = self.LRiterator(self)
        return iterator
    def add_element(self,element):
        """ add element to lattice """
        if len(self.seq) == 0:
            s0 = 0.
            element.prev = None
        else:
            previous = self.seq[-1]  
            previous.next = element   # forward link
            element.prev  = previous  # back link
            s0 = element.prev.position[2]  # end of previous
        self.seq.append(element)
        l = element.length
        si = s0
        sf = si+l
        sm = (sf+si)/2.
        position = (si,sm,sf)
        element.position = position 
        self.length = sf
    def toString(self):
        # TODO needs improvement
        """ log lattice layout to string (could be even better?) """
        mcell = ELM.I(label='')   #  chain matrices
        # for element in self.seq:
        for element in iter(self):
        # for element in iter(self):
            DEBUG_OFF('{:10s}({:d})\tlength={:.3f}\tfrom-to: {:.3f} - {:.3f}'.format(element.label,id(element),element.length,element.position[0],element.position[2]))
            # ACHTUNG: Reihenfolge im Produkt ist wichtig!
            mcell = element * mcell   
        mcell.section = '<= full lattice map'
        return mcell.prmatrix()
    def stats(self,soll_track):
        """ gather lattice statistics """
        cav_counter = 0
        q_counter   = 0
        ttfm = +1.e+50
        ttfx = +1.e-50
        tk_i = soll_track.getpoints()[0]()[6]
        tk_f = soll_track.getpoints()[-1]()[6]
        """ loop over all nodes in lattice """
        for element in iter(self):
            if isinstance(element,(ELM.QF,ELM.QD)):
                q_counter += 1
            if isinstance(element,(ELM.RFG,ELM.RFC)):
                cav_counter += 1
                ttfm = min(element.ttf,ttfm)
                ttfx = max(element.ttf,ttfx)
        if q_counter == 0:
            SUMMARY['nbof quadrupoles*'] = '0 (no thick quads?)'
        else:
            SUMMARY['nbof quadrupoles*'] = q_counter
        SUMMARY['nbof cavities*']        = cav_counter
        SUMMARY['(ttf)min,(ttf)max*']    = (ttfm,ttfx)
        SUMMARY['(energy)i,(energy)f [MeV]']  = (tk_i,tk_f)
    def cell(self,closed=True):
        """
        Construct the full accelerator lattice-cell matrix and extract standard quantities:
            full cell: mcell
            stability?
            betatron tunes: mux, muy
            det(M)
            check symplecticity
            twiss prameters beta, alpha, gamma for periodic lattices
        """
        mcell = ELM.I(label=' <==')   #  chain matrices for full cell
        """ loop over all lattice nodes """
        for count,element in enumerate(iter(self)):
            if count == 0:
                mcell = element
            else:
                # Achtung: Reihenfolge !
                mcell = element * mcell

        ## Stabilität ?
        unstable = False
        stab = fabs(mcell.tracex())
        PARAMS['traceX'] = stab
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
        PARAMS['traceY'] = stab
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
        self.accel = mcell    # the full cell: isinstance(self.accel,Lattice)==True
        # if verbose:
        printv(0,'Full Accelerator Matrix (f)<==(i)')
        printv(0,self.accel.prmatrix())
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
        emitx = PARAMS['emitx_i']
        emity = PARAMS['emity_i']

        ## Startwerte für twiss-functions aus cell matrix (not beta_matrix!)
        if closed:
            if not unstable:
                cell_matrix = self.accel.matrix
                m11  = cell_matrix[XKOO,XKOO];   m12  = cell_matrix[XKOO,XPKOO]
                m21  = cell_matrix[XPKOO,XKOO];  m22  = cell_matrix[XPKOO,XPKOO]
                n11  = cell_matrix[YKOO,YKOO];   n12  = cell_matrix[YKOO,YPKOO]
                n21  = cell_matrix[YPKOO,YKOO];  n22  = cell_matrix[YPKOO,YPKOO]
                ## betax,alphax,gammax from transfer matrix
                #  [m11,m12] = [cos(mu) + alpha * sin(mu),       beta * sin(mu)     ]
                #  [m21,m22] = [     -gamma * sin(mu)    , cos(mu) - alpha * sin(mu)]
                #  und: beta * gamma - alpha^2 = 1
                sin2mu = -(((m11-m22)**2)/4.+m12*m21)
                if sin2mu < 0.: 
                    print('STOP: unstable in X: sin(mux)^2={:4.4f}'.format(sin2mu))
                    sys,exit(1)
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
                if sin2mu < 0.: 
                    print('STOP: unstable in Y: sin(muy)^2={:4.4f}'.format(sin2mu))
                    sys.exit(1)
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
                v_beta_a = NP.array([bax,alx,gmx,bay,aly,gmy,1.,0.,1.])
                m_cell_beta = self.accel.beta_matrix()
                v_beta_e = NP.dot(m_cell_beta,v_beta_a)
                # if verbose:
                printv(1,'Probe: {TW(f)} == {BetaMatrix}x{TW(i)}?')
                diffa_e = v_beta_a - v_beta_e
                for i in range(6):
                    if fabs(diffa_e[i]) < 1.e-9: diffa_e[i] = 0.
                printv(1,'TW(i)-TW(f) (should be [0,...,0]):\n',diffa_e)
                ## transversale twiss parameter fuer periodisches lattice
                twx = Twiss(bax,alx,emitx)
                twy = Twiss(bay,aly,emity)
                PARAMS['twiss_x_i'] = twx
                PARAMS['twiss_y_i'] = twy
                SUMMARY['(sigx )i*   [mm]'] =  PARAMS['twiss_x_i'].sigmaH()*1.e3
                SUMMARY["(sigx')i* [mrad]"] =  PARAMS['twiss_x_i'].sigmaV()*1.e3
                SUMMARY['(sigy )i*   [mm]'] =  PARAMS['twiss_y_i'].sigmaH()*1.e3
                SUMMARY["(sigy')i* [mrad]"] =  PARAMS['twiss_y_i'].sigmaV()*1.e3
            else:
                print('STOP: unstable lattice!')
                sys.exit(1)
        else:
            ## transversale twiss parameter fuer transfer lines
            # alfa, beta und emittance definieren den beam @ entrance, 
            # NOTE: transfer lattices need not to be stable!
            bax,alx,gmx,epsx = PARAMS['twiss_x_i']()
            bay,aly,gmy,epsy = PARAMS['twiss_y_i']()
        printv(0,'using @ entrance: [beta,  alfa,  gamma]-X    [beta,   alfa,   gamma]-Y')
        printv(0,'                  [{:.3f}, {:.3f}, {:.3f}]-X    [{:.3f},  {:.3f},  {:.3f}]-Y'.format(bax,alx,gmx,bay,aly,gmy))
        return(self)
    def report(self):
        # TODO needs more work
        """ report lattice layout (may not work!) """
        raise RuntimeWarning('Lattice.report() not ready')
        reprt = ''
        header = ''
        row = ''
        # for count,element in enumerate(reversed(self.seq)):
        for count,element in enumerate(iter(self)):
            name = element.label
            len = element.length
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
    def toggle_iteration(self):
        """ toggle l->R or L<-R sweep through linked list of elements"""
        if self.iteration == "LR":
            self.iteration = "RL"
        elif self.iteration == "RL":
            self.iteration = "LR"
    def concat(self,lattice):
        """Concatenate two Lattice pieces (self+lattice)"""
        # TODO check correctness again
        for element in iter(lattice):
            element = copy(element)
            self.add_element(element)
    def twiss_envelopes(self,steps=1):
        """ Calulate envelopes from initial twiss-vector with beta-matrices """
        twfun = Functions(('s','bx','ax','gx','by','ay','gy','bz','az','gz'))
        bx,ax,gx,epsx = PARAMS['twiss_x_i']()
        by,ay,gy,epsy = PARAMS['twiss_y_i']()
        bz,az,gz,epsz = PARAMS['twiss_z_i']()
        twiss_vector0 = NP.array([bx,ax,gx,by,ay,gy,bz,az,gz])   # initial
        B_matrix = NP.eye(9)                            # cumulated beta-matrix
        """ loop over all nodes in the lattice """
        for node in iter(self):
            slices = node.make_slices(anz = steps)
            means = []
            s,sm,sf = node.position
            for slice in slices:
                beta_matrix = slice.beta_matrix()
                B_matrix    = NP.dot(beta_matrix,B_matrix)
                twiss_vector= NP.dot(B_matrix,twiss_vector0)     # track twiss-vector
                s += slice.length
                twfun.append(s,tuple(twiss_vector))
                node['twiss'] = tuple(twiss_vector)
                bx    = twiss_vector[Ktw.bx]; ax = twiss_vector[Ktw.ax]; gx = twiss_vector[Ktw.gx]
                by    = twiss_vector[Ktw.by]; ay = twiss_vector[Ktw.ay]; gy = twiss_vector[Ktw.gy]
                sigxy = (*sigmas(ax,bx,epsx),*sigmas(ay,by,epsy))
                means.append(sigxy)
                
            # means = NP.array(means)
            means = NP.mean(means,axis=0)
            # each node has its tuple of average sigmas
            node['sigxy'] = tuple(means)
            # aperture check
            self.aperture_check(node,twiss=True)
        return twfun
    def aperture_check(self,node,twiss=True):
        """ check sigmas against apertures """
        fcnt = 'twiss_envelopes()' if twiss else 'sigma_envelopes()'
        s,sm,sf = node.position
        if FLAGS['useaper']:
            nbsigma = PARAMS['nbsigma']
            if node.aperture != None:
                aperture = node.aperture
                sigx, sigxp, sigy, sigyp = node['sigxy']
                if PARAMS['warnmx']:
                    if(aperture < nbsigma*sigx or aperture < nbsigma*sigy):
                        # warnings.showwarning('{} sigma aperture hit @ s={:.1f} [m]'.format(nbsigma,sm),UserWarning,'lattice.py',fcnt)
                        print(warnings.formatwarning('{}: {} sigma aperture hit @ s={:.1f} [m]'.format(fcnt,nbsigma,sm),UserWarning,'','')[3:-1])
                        PARAMS['warnmx'] -= 1
                        if PARAMS['warnmx'] == 0: print('skipping more warnings ...')
    def sigmas(self,steps = 10):
        """ dispatch to different envelope functions """
        #TODO: analytical sigmas for 'dyn' mapping as best estimates ?
        def envelopes(function, steps = 10):
            """ calc. envelopes using function """
            twfunc     = function(steps = steps)
            b,a,g,epsx = PARAMS['twiss_x_i']()
            b,a,g,epsy = PARAMS['twiss_y_i']()
            sigma_fun  = Functions(('s','sigmax','sigmay'))
            for i in range(twfunc.nbpoints):
                val,dummy = twfunc[i]
                (s,bx,ax,gx,by,ay,gy,bz,az,gz) = val
                val=(sqrt(bx*epsx),sqrt(by*epsy))
                sigma_fun.append(s,val)
            return sigma_fun

        if PARAMS['mapping'] == 'dyn':
            mess = 'CALCULATE TWISS ENVELOPES WITH T3D CAVITIES ("dyn")'
            function = self.twiss_envelopes # use beta-matrix            
        elif FLAGS['sigma']:
            mess = 'CALCULATE SIGMA ENVELOPES'
            function = self.sigma_envelopes # use sigma-matrix
        elif not FLAGS['sigma']:
            mess = 'CALCULATE TWISS ENVELOPES'
            function = self.twiss_envelopes # use beta-matrix

        if not FLAGS['KVout']: 
            print(mess)
            sigma_fun = envelopes(function, steps = steps)
        return sigma_fun
    def sigma_envelopes(self, steps = 1):
        """ Envelopes and twiss-functions from sigma-matrix method a.k.a rms-envelopes """
        # initials
        bx,ax,gx,epsx = PARAMS['twiss_x_i']()
        by,ay,gy,epsy = PARAMS['twiss_y_i']()
        bz,az,gz,epsz = PARAMS['twiss_z_i']()
        twiss_vector0     = NP.array([bx,ax,gx,by,ay,gy,bz,az,gz])  # twiss vector IN lattice
        sg0      = Sigma(twiss_vector0,epsx,epsy,epsz)              # sigma object IN lattice
        # sigma envelopes as function of distance s
        sigma_fun = Functions(('s','bx','ax','gax','by','ay','gy','bz','az','gz'))
        for node in iter(self): # loop nodes
            # sigma-matrices for a single node
            sigmas = node.sigma_beam(steps = steps, sg = sg0)   # note: sets node['sigxy']
            # prep plot list of ftn's
            for sg,s in sigmas:
                v = sg.twiss()      # twiss from Sigma object
                flist = v.tolist()
                sigma_fun.append(s,tuple(flist))
            sg0 = sg          # loop back nodes

            # aperture check
            self.aperture_check(node,twiss=False)
        return sigma_fun
    def dispersion(self,steps=10,closed=True):
        """ track the dispersion function """
        traj = []
        # column vector with MDIM rows, 1 column
        v_0 = NP.array([0.,0.,0.,0.,0.,1.,0.,0.,0.,0.])
        if closed == True:
            m_cell = self.accel
            m11 = m_cell.matrix[0,0]
            m15 = m_cell.matrix[0,5]
            d0  =  m15/(1.-m11)     # from H.Wiedemann (6.79) pp.206
            v_0[0] = d0
        s = 0.0
        traj = [(s,v_0[0],v_0[1])]
        # for element in self.seq:
        for element in iter(self):
            slices = element.make_slices(anz = steps)
            for i_element in slices:
                m_beta = i_element.matrix
                v_0 = NP.dot(m_beta,v_0)
                s += i_element.length
                d  = v_0[0]
                dp = v_0[1]
                traj.append((s,d,dp))
        return traj
    def lattice_plot_functions(self):
        """
        generate the functions to plot the lattice
        viseo:    shows elements
        aperture: shows the physical apertures of elements
        """
        fun = Functions(('s','viseo'))
        ape = Functions(('s','aperture'))
        # for element in self.seq:
        for element in iter(self):
            # DEBUG((element.__class__,element['viseo'],element.position))
            pos   = element.position
            # element plot
            viseo = element['viseo']
            si, sm, sf = pos
            fun.append(si,(0,))
            fun.append(si,(viseo,))
            fun.append(sf,(viseo,))
            fun.append(sf,(0,))

            # aperture plot
            aperture = None
            if hasattr(element, 'aperture'): 
                aperture = element.aperture
            if element.__class__.__name__ == 'D': continue
            if aperture == None: continue
            ape.append(sm,(aperture,))
        return fun,ape
    def cs_traj(self,steps=10):
        """ track cos- & sin-trajectories """
        def SollTest_ON(arg):  # set all 0. to simulate Sollteilchen
            return 0.
        def SollTest_OFF(arg):
            return arg
        soll_test = SollTest_OFF

        def DEBUG_TRACKs(elm,c,s):
            DEBUG_ON()
            print('{} cosine {} sine {}'.format(elm.type,c,s))
            pass
        # function body --------------- function body --------------- function body --------------- 
        # function body --------------- function body --------------- function body --------------- 
        # function body --------------- function body --------------- function body --------------- 
        print('CALCULATE C+S TRAJECTORIES')
        tkin = PARAMS['sollteilchen'].tkin
        
        """ injektion parameters """
        if True:
            # 2 point on the ellipse y1 & y4: intersections
            x1,x1p = soll_test(PARAMS['twiss_x_i'].y1())
            y1,y1p = soll_test(PARAMS['twiss_y_i'].y1())
            x4,x4p = soll_test(PARAMS['twiss_x_i'].y4())
            y4,y4p = soll_test(PARAMS['twiss_y_i'].y4())
        else:
            # 2 point on the ellipse y2 & y3: maximum values
            x1,x1p = soll_test(PARAMS['twiss_x_i'].y2())
            y1,y1p = soll_test(PARAMS['twiss_y_i'].y2())
            x4,x4p = soll_test(PARAMS['twiss_x_i'].y3())
            y4,y4p = soll_test(PARAMS['twiss_y_i'].y3())
        if FLAGS['dWf']:
            sigmaz_i    = soll_test(PARAMS['z0'])      # z0[m] from waccept
            Dp2p_i      = soll_test(PARAMS['Dp2p0'])   # dp/p0 from waccept
        else:
            sigmaz_i = 0.
            Dp2p_i   = 0.
        z1,z1p = soll_test((sigmaz_i, 0.))    # S
        z4,z4p = soll_test((0., Dp2p_i))      # C
        # MDIMxMDIM tracking used here
        s   = 0.
        c_0 = NP.array([x1, x1p, y1, y1p, z4, z4p, tkin,1,0,1])  # C
        s_0 = NP.array([x4, x4p, y4, y4p, z1, z1p, tkin,1,0,1])  # S
        # function names
        c_fun = Functions(('s','cx','cxp','cy','cyp','cz','cdp'))
        s_fun = Functions(('s','sx','sxp','sy','syp','sz','sdp'))

        """ loop through lattice """
        for element in iter(self):
            # objprnt(particle,text='cs_traj: '+element.label)   # DEBUG
            slices = element.make_slices(anz=steps)
            try:
                for i_element in slices:
                    # if i_element.type == 'QFth' : print(i_element.type,i_element.matrix)
                    # if i_element.type == 'QDth' : print(i_element.type,i_element.matrix)
                    s += i_element.length
                    ## COSine_like
                    c_0 = i_element.map(c_0)   # map!!!
                    cx  = c_0[XKOO]
                    cxp = c_0[XPKOO]
                    cy  = c_0[YKOO]
                    cyp = c_0[YPKOO]
                    # cz  = -c_0[ZKOO]*360./(beta*lamb)            # sigmaz_i --> dPhi [deg]
                    # cdw = c_0[ZPKOO]*(gamma+1.)/gamma*100.       # dp/p --> dW/W [%]
                    cz  = c_0[ZKOO]*1.e3      # z [mm]
                    cdp = c_0[ZPKOO]*100.     # dp/p [%]
                    c_fun.append(s,(cx,cxp,cy,cyp,cz,cdp))
                    ## SINus_like
                    s_0 = i_element.map(s_0)   # map!!!
                    sx  = s_0[XKOO]
                    sxp = s_0[XPKOO]
                    sy  = s_0[YKOO]
                    syp = s_0[YPKOO]
                    # sz  = -s_0[ZKOO]*360./(beta*lamb)
                    # sdw = s_0[ZPKOO]*(gamma+1.)/gamma*100.
                    sz  = -s_0[ZKOO]*1.e3
                    sdp = s_0[ZPKOO]*100.
                    s_fun.append(s,(sx,sxp,sy,syp,sz,sdp))
                    if 0: DEBUG_TRACKs(i_element,(cx,cxp,cy,cyp,cz,cdp),(sx,sxp,sy,syp,sz,sdp))
            except ValueError as ex:
                print('STOP: C+S TRAJECTORIES at s = {:6.2f} [m]!'.format(s))
                #TODO  raise ex
                break
        return (c_fun,s_fun)
    def symplecticity(self):
        """ test symplecticity """
        s = NP.array([
                    [ 0.,1., 0.,0., 0.,0.,0.,0.,0.,0.],    #x
                    [-1.,0., 0.,0., 0.,0.,0.,0.,0.,0.],    #x'
                    [ 0.,0., 0.,1., 0.,0.,0.,0.,0.,0.],    #y
                    [ 0.,0.,-1.,0., 0.,0.,0.,0.,0.,0.],    #y'
                    [ 0.,0., 0.,0., 0.,1.,0.,0.,0.,0.],    #z
                    [ 0.,0., 0.,0.,-1.,0.,0.,0.,0.,0.],    #z'
                    [ 0.,0., 0.,0., 0.,0.,1.,0.,0.,0.],    #T
                    [ 0.,0., 0.,0., 0.,0.,0.,1.,0.,0.],    #1
                    [ 0.,0., 0.,0., 0.,0.,0.,0.,1.,0.],    #S
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
    def show_linkage(self):
        """ Show left-right links of doubly linked element list of the lattice. Iterate in both directions """
        print("@@@@@@@@@@ iteration {:s} @@@@@@@@@@".format(self.iteration))
        for next in iter(self):
            print('{:38s}\t{:38s}'.format(repr(next),repr(next.next)))
            next = next.next
        self.toggle_iteration()
        print("@@@@@@@@@@ iteration {:s} @@@@@@@@@@".format(self.iteration))
        for next in iter(self):
            print('{:38s}\t{:38s}'.format(repr(next),repr(next.prev)))
            next = next.next
        self.toggle_iteration()
    @property
    def first_gap(self):
        """ return the 1st RF gap"""
        node = None
        # for elm in self.seq:
        for elm in iter(self):
            if isinstance(elm,(ELM.RFG,ELM.RFC,ELM.GAP)):
                node = elm
                break
        return node
## utilities
def make_wille():
    """
    Wille's test lattice
    """
    print("K.Wille's Beispiel auf pp.113 Formel (3.200)")
    kqf = wille()['k_quad_f']
    lqf = wille()['length_quad_f']
    kqd = wille()['k_quad_d']
    lqd = wille()['length_quad_d']
    rhob = wille()['bending_radius']
    lb = wille()['dipole_length']
    ld = wille()['drift_length']
    # elements
    mqf1 = ELM.QF(k0=kqf,length=lqf,label='QF1')
    mqf2 = ELM.QF(k0=kqf,length=lqf,label='QF2')
    mqd1 = ELM.QD(k0=kqd,length=lqd,label='QD1')
    md1  = ELM.D(length=ld)
    md2  = ELM.D(length=ld)
    md3  = ELM.D(length=ld)
    md4  = ELM.D(length=ld)
    mbr1  = ELM.RD(radius=rhob,length=lb)
    mbr2  = ELM.RD(radius=rhob,length=lb)
    # lattice
    lattice = Lattice()
    lattice.add_element(mqf1)
    lattice.add_element(md1)
    lattice.add_element(mbr1)
    lattice.add_element(md2)
    lattice.add_element(mqd1)
    lattice.add_element(md3)
    lattice.add_element(mbr2)
    lattice.add_element(md4)
    lattice.add_element(mqf2)
    DEB.get('OFF')('lattice: {}'.format(lattice.toString()))
    top = Lattice()
    top.concat(lattice)
    top.concat(lattice)
    top.concat(lattice)
    return top
def test1():
    from matplotlib.pyplot import plot,show,legend
    print('-------------------------------------Test1--')
    lattice = make_wille()
    # for element in lattice.seq: DEBUG('test1: ',' [si,sm,sf,] elm= [{1}] {0}'.format(repr(element),''.join('{:5.3f},'.format(el) for el in element.position)))
    # cell boundaries
    full_cell = lattice.cell(closed=True)
    betax,a,g,e = PARAMS['twiss_x_i']()
    betay,a,g,e = PARAMS['twiss_y_i']()
    PARAMS['twiss_z_i'] = Twiss(1.,0.,1.)
    lattice.symplecticity()
    # twiss functions
    beta_fun = lattice.twiss_envelopes(steps=5)
    # cl,sl = lattice.cs_traj(steps=100)sK
    disp = lattice.dispersion(steps=100,closed=True)
    # plots
    s  = [beta_fun(i,'s')  for i in range(beta_fun.nbpoints)]
    xs = [beta_fun(i,'bx') for i in range(beta_fun.nbpoints)]
    ys = [beta_fun(i,'by') for i in range(beta_fun.nbpoints)]
    sd = [x[0] for x in disp]             # abszisse s
    ds = [x[1] for x in disp]             # dispersion(s)
    #-------------------- lattice viseo
    lat_plot, ape_plot = lattice.lattice_plot_functions()
    vis_abszisse = [lat_plot(i,'s')              for i in range(lat_plot.nbpoints)]
    vis_ordinate = [lat_plot(i,'viseo')          for i in range(lat_plot.nbpoints)]
    vzero        = [0.                           for i in range(lat_plot.nbpoints)] # zero line

    plot(s,xs,label='betax')
    plot(s,ys,label='betay')
    plot(sd,ds,label='disp')
    plot(vis_abszisse,vis_ordinate,label='',color='black')
    plot(vis_abszisse,vzero,color='black')
    legend(loc='upper left')
    show()
def test2():
    print('-------------------------------------Test2--')
    print('lattice tags test ...')
    lattice = Lattice()
    print(type(lattice))
    lattice.name = 'NAME'
    lattice.label = 'LABEL'
    lattice.section = 'SECTION'
    print(lattice.__dict__)
def test3():
    print('-------------------------------------Test3--')
    lattice = make_wille()
    lattice.show_linkage()
if __name__ == '__main__':
    test1()
    test2()
    test3()
