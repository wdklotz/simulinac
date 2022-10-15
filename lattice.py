#!/Users/klotz/SIMULINAC_env/bin/python
# -*- coding: utf-8 -*-
__version__='v10.22.3'
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
import unittest
from math import sqrt,fabs,acos,degrees
from numpy import linalg as LA
from copy import copy
# from termcolor import colored
# from sty import fg,bg,ef,rs

import elements as ELM
from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO
from setutil import PARAMS,FLAGS,SUMMARY,print_verbose,sigmas, objprnt, Ktw, Ktp
from setutil import Twiss, Functions, Particle, Proton, colors, MDIM, DEBUG_ON, DEBUG_OFF
from sigma import Sigma, sig_map

class Lattice(object):
    """ The Lattice object is a list of elements: ELM.<element> in self.seq  ?? """
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
    def __init__(self, injection_energy=50.,descriptor=""):
        self.seq              = []       # list of _Node objects z.B. [D,QD,GAP,QF....]
        self.iteration        = "LR"  # default: iterating lattice left-right
        self.length           = 0.
        self.acc_node         = None
        self.injection_energy = injection_energy
        self.descriptor       = descriptor
    def __iter__(self):
        """ iterator using the linked list of element """
        if self.iteration == "RL":
            iterator = self.RLiterator(self)
        elif self.iteration == "LR":
            iterator = self.LRiterator(self)
        return iterator
    def add_node(self,node):
        """ 
        Add node to end of lattice. lattice orientation from left to right.
        Track the reference particle to get the energy kicks and length advances using t3d matrices,.
        Link the lattice nodes in a doubly linked list.
        Calculate the node positions.
         """
        """ the 1st node """
        if len(self.seq) == 0:
            """ the 1st node """
            tk_injection = self.injection_energy
            ref_track    = NP.array([0.,0.,0.,0.,0.,0.,tk_injection,1.,0.,1.])
            node_adj     = node.adjust_energy(tk_injection)    # energy ADJUST
            # ref_track    = NP.dot(node.matrix,ref_track)     # track @ out of 1st node
            ref_track_m  = node_adj.map(ref_track)
            ref_particle = Proton(ref_track_m[EKOO])           # ref_particle @ out of 1st node
            node_adj.prev = node.next = None
            sf = node_adj.length
            sf = ref_track_m[SKOO]
            node_adj.position     = (0.,sf/2.,sf)
            node_adj.ref_track    = ref_track_m
            node_adj.ref_particle = ref_particle
        else:
            """ all nodes after 1st """
            prev         = self.seq[-1]
            ref_track    = prev.ref_track
            si           = ref_track[SKOO]
            tkin         = ref_track[EKOO]
            node_adj     = node.adjust_energy(tkin)
            # ref_track = NP.dot(node.matrix,ref_track)
            ref_track_m  = node_adj.map(ref_track)
            ref_particle = Proton(ref_track_m[EKOO])
            prev.next = node_adj
            node_adj.prev = prev
            sf        = ref_track_m[SKOO]
            node_adj.position     = (si,(si+sf)/2,sf)
            node_adj.ref_track    = ref_track_m
            node_adj.ref_particle = ref_particle
        self.length = sf    # lattice length
        self.seq.append(node_adj)
    def toString(self):
       # TODO needs improvement
        """ log lattice layout to string (could be even better?) """
        mcell = ELM.I(label='')   #  chain matrices
        # for element in self.seq:
        for element in iter(self):
        # for element in iter(self):
            DEBUG_OFF('{:10s}({:d}) length={:.3f} from-to: {:.3f} - {:.3f}'.format(element.label,id(element),element.length,element.position[0],element.position[2]))
            # ACHTUNG: Reihenfolge im Produkt ist wichtig!
            mcell = element * mcell   
        mcell.section = '<= full lattice map'
        return mcell.prmatrix()
    def stats(self):
        """ gather lattice statistics """
        cavity_counter = 0
        quad_counter   = 0
        ttfmin = +1.e+50
        ttfmax = +1.e-50
        # tk_i = soll_track.getpoints()[0]()[6]
        # tk_f = soll_track.getpoints()[-1]()[6]
        tk_inject  = self.injection_energy
        tk_final   = self.seq[-1].ref_track[EKOO]
        """ loop over all nodes in lattice """
        for element in iter(self):
            if isinstance(element,(ELM.QF,ELM.QD)):
                quad_counter += 1
            if isinstance(element,(ELM.RFG,ELM.RFC)):
                cavity_counter += 1
                ttfmin = min(element.ttf,ttfmin)
                ttfmax = max(element.ttf,ttfmax)
        if quad_counter == 0:
            SUMMARY['nbof quadrupoles*'] = '0 (no quads?)'
        else:
            SUMMARY['nbof quadrupoles*']      = quad_counter
        SUMMARY['nbof cavities*']             = cavity_counter
        SUMMARY['(ttf)min,(ttf)max*']         = (ttfmin,ttfmax)
        SUMMARY['(energy)i,(energy)f [MeV]']  = (tk_inject,tk_final)
        SUMMARY['lattice length [m]']         = self.length
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
        # mcell = ELM.I(label=' <==')   #  chain matrices for full cell
        mcell = None
        """ loop over all lattice nodes """
        for count,element in enumerate(iter(self)):
            if count == 0:
                mcell = element
            else:
                # Achtung: Reihenfolge !
                mcell = element * mcell    
                # mcell = mcell *element   #???????

        ## Stabilität ?
        unstable = False
        stab = fabs(mcell.tracex())
        # if verbose:
        print_verbose(1,'stability X? ',stab)
        if stab >= 2.0:
            # if verbose:
            print_verbose(1,'unstable Lattice in x-plane\n')
            unstable = True
        else:
            cos_mux = 0.5 * stab
            mux = degrees(acos(cos_mux))

        stab = fabs(mcell.tracey())
        # if verbose:
        print_verbose(1,'stability Y? ',stab)
        if stab >= 2.0:
            # if verbose:
            print_verbose(1,'unstable Lattice in y-plane\n')
            unstable = True
        else:
            cos_muy = 0.5 * stab
            muy = degrees(acos(cos_muy))
        if not unstable:
            # if verbose:
            print_verbose(1,'\nphase_advance: X[deg]={:3f} Y[deg]={:.3f}\n'.format(mux,muy))
        ## full accelerator
        self.acc_node = mcell    # the full acc: isinstance(self.acc_node,Lattice)==True
        # if verbose:
        print_verbose(0,'Full Accelerator Matrix (f)<==(i)')
        print_verbose(0,self.acc_node.prmatrix())
        det = LA.det(self.acc_node.matrix)
        # if verbose:
        print_verbose(2,'det|full-cell|={:.5f}\n'.format(det))
        ## Determinate M-I == 0 ?
        beta_matrix = mcell.beta_matrix()
        for i in range(5):
            beta_matrix[i,i] = beta_matrix[i,i]-1.0
        det = LA.det(beta_matrix)
        # if verbose:
        print_verbose(2,'det|Mbeta - I|={:.5f}\n'.format(det))
        ## symplectic?
        s = self.symplecticity()
        # if verbose:
        print_verbose(2,'symplectic (+1,-1,+1,-1,+1,-1)?')
        print_verbose(2,'[{:4>+.2f}, {:4>+.2f}, {:4>+.2f}, {:4>+.2f}, {:4>+.2f}, {:4>+.2f}]\n'.format(s[0],s[1],s[2],s[3],s[4],s[5]))
        ## Vorgabe emittance @ entrance
        emitx = PARAMS['emitx_i']
        emity = PARAMS['emity_i']

        ## Startwerte für twiss-functions aus cell matrix (not beta_matrix!)
        if closed:
            if not unstable:
                cell_matrix = self.acc_node.matrix
                m11  = cell_matrix[XKOO,XKOO];   m12  = cell_matrix[XKOO,XPKOO]
                m21  = cell_matrix[XPKOO,XKOO];  m22  = cell_matrix[XPKOO,XPKOO]
                n11  = cell_matrix[YKOO,YKOO];   n12  = cell_matrix[YKOO,YPKOO]
                n21  = cell_matrix[YPKOO,YKOO];  n22  = cell_matrix[YPKOO,YPKOO]
                # Formeln von Wille pp 107
                wurzx = 2.-m11**2-2.*m12*m21-m22**2
                if wurzx > 0.:
                    bax = fabs(2.*m12/sqrt(wurzx))     #  muss Absolutwert sein! (Fehler bei Wille)
                    alx = (m11-m22)/(2.*m12)*bax
                    gmx = (1.+alx**2)/bax
                print('betax {:4.4f} alfax {:4.4f} gammax {:4.4f}'.format(bax,alx,gmx))
                wurzy = 2.-n11**2-2.*n12*n21-n22**2
                if wurzy > 0.:
                    bay = fabs(2.*n12/sqrt(wurzy))    #  muss Absolutwert sein! (Fehler bei Wille)
                    aly = (n11-n22)/(2.*n12)*bay
                    gmy = (1.+aly**2)/bay
                print('betay {:4.4f} alfay {:4.4f} gammay {:4.4f}'.format(bay,aly,gmy))
                ## Probe: twiss-functions durch ganze Zelle mit beta-matrix (nur sinnvoll fuer period. Struktur!)
                v_beta_a = NP.array([bax,alx,gmx,bay,aly,gmy,1.,0.,1.])
                m_cell_beta = self.acc_node.beta_matrix()
                v_beta_e = NP.dot(m_cell_beta,v_beta_a)
                # if verbose:
                print_verbose(1,'Probe: {TW(f)} == {BetaMatrix}x{TW(i)}?')
                diffa_e = v_beta_a - v_beta_e
                for i in range(6):
                    if fabs(diffa_e[i]) < 1.e-9: diffa_e[i] = 0.
                print_verbose(1,'TW(i)-TW(f) (should be [0,...,0]):\n',diffa_e[0:6])
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
            # alfa, beta und emittance definieren den Beam @ entrance, 
            bax,alx,gmx,epsx = PARAMS['twiss_x_i']()
            bay,aly,gmy,epsy = PARAMS['twiss_y_i']()
        print_verbose(0,'using @ entrance: [beta,  alfa,  gamma]-X    [beta,   alfa,   gamma]-Y')
        print_verbose(0,'                  [{:.3f}, {:.3f}, {:.3f}]-X    [{:.3f},  {:.3f},  {:.3f}]-Y'.format(bax,alx,gmx,bay,aly,gmy))
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
    def twiss_funcs(self,steps=1):
        """ Calulate twiss functions and evelope functions from initial twiss-vector with beta-matrices or sigma-matrix"""
        sFLAG   = FLAGS['sigma']
        nlFLAG  = FLAGS['non_linear_mapping']
        t3dFLAG = True
        
        mess = ""
        if nlFLAG:
            mess = colors.RED+'WARN: Lattice has RF-gaps with non-linear mapping. ENVELOPES are calulated using T3D\'s RF-gaps (NT=10) instead.\n'+colors.ENDC
            mess += 'sigma ENVELOPES from TWISS paramweters'
        elif not nlFLAG:
            if sFLAG:
                mess += 'sigma ENVELOPES from SIGMA-matrix formalism'
                t3dFLAG = False
            else:
                mess += 'sigma ENVELOPES from TWISS paramweters'
        if not FLAGS['KVout']: 
            print(mess)
        
        # INITIAL twiss @ entrance
        function_tbl = []
        bx,ax,gx,epsx = PARAMS['twiss_x_i']()
        by,ay,gy,epsy = PARAMS['twiss_y_i']()
        bz,az,gz,epsz = PARAMS['twiss_z_i']()
        node0 = self.seq[0]      # 1st node
        si,sm,sf      = node0.position
        twiss_vector0 = NP.array([bx,ax,gx,by,ay,gy,bz,az,gz])
        sigx,sigxp,sigy,sigyp = (*sigmas(ax,bx,epsx),*sigmas(ay,by,epsy))
        node0.twiss  = tuple(twiss_vector0)
        node0.sigxy  = (sigx,sigxp,sigy,sigyp)
        function_row = (si,bx,ax,gx,by,ay,gy,bz,az,gz,sigx,sigxp,sigy,sigyp)
        function_tbl.append(function_row)
                
        """ loop over all nodes in the lattice to get sliced function values"""
        B_matrix = NP.eye(9,9)                           # cumulated beta-matrix
        R_matrix = NP.eye(MDIM,MDIM)                     # cumulated R-matrix
        Sig      = Sigma(twiss_vector0,epsx,epsy,epsz)   # cumulated Sigma object
        for node in iter(self):
            si,sm,sf = node.position
            slices = node.make_slices(anz=steps)
            s = si
            means = []
            for slice in slices:
                slice_beta_mx = slice.beta_matrix()
                slice_r_mx    = slice.matrix
                B_matrix      = NP.dot(slice_beta_mx,B_matrix)
                R_matrix      = NP.dot(slice_r_mx,R_matrix)
                Sig           = sig_map(Sig,slice)
                
                if t3dFLAG:
                    twiss_vector = NP.dot(B_matrix,twiss_vector0)     # track twiss-vector
                    sigx,sigxp,sigy,sigyp = (*sigmas(ax,bx,epsx),*sigmas(ay,by,epsy))
                elif not t3dFLAG:
                    twiss_vector = Sig.sig_twiss_vec_get()
                    sigma_vector = Sig.sig_sigma_vec_get()
                    sigx  = sigma_vector[0]; sigxp = sigma_vector[1]              
                    sigy  = sigma_vector[2]; sigyp = sigma_vector[3]              
                    # sigz  = sigma_vector[4]; sigzp = sigma_vector[5]              
                
                bx = twiss_vector[Ktw.bx]; ax = twiss_vector[Ktw.ax]; gx = twiss_vector[Ktw.gx]
                by = twiss_vector[Ktw.by]; ay = twiss_vector[Ktw.ay]; gy = twiss_vector[Ktw.gy]
                s += slice.length
                function_row = (s,bx,ax,gx,by,ay,gy,bz,az,gz,sigx,sigxp,sigy,sigyp)
                function_tbl.append(function_row)
                means.append((sigx,sigxp,sigy,sigyp))
            means = NP.mean(means,axis=0)
            node.twiss = tuple(twiss_vector)    # each noe has twiss
            node.sigxy = tuple(means)           # each node has sigxy
            # aperture check
            self.aperture_check(node,twiss=t3dFLAG)
            
        twissfun = Functions(('s','bx','ax','gx','by','ay','gy','bz','az','gz','sigx','sigxp','sigy','sigyp'))
        for row in function_tbl:
            abscisse  = row[0]
            ordinaten = row[1:]
            twissfun.append(abscisse,ordinaten)
        return twissfun
    def aperture_check(self,node,twiss=True):
        """ check sigmas against apertures """
        fcnt = 'twiss envelopes' if twiss else 'sigma envelopes'
        s,sm,sf = node.position
        if FLAGS['useaper']:
            nbsigma = PARAMS['nbsigma']
            if node.aperture != None:
                aperture = node.aperture
                sigx, sigxp, sigy, sigyp = node.sigxy
                if PARAMS['warnmx']:
                    if(aperture < nbsigma*sigx or aperture < nbsigma*sigy):
                        # warnings.showwarning('{} sigma aperture hit @ s={:.1f} [m]'.format(nbsigma,sm),UserWarning,'lattice.py',fcnt)
                        print(warnings.formatwarning(colors.RED+'{}: {} sigma aperture hit @ s={:.1f} [m]'.format(fcnt,nbsigma,sm)+colors.ENDC,UserWarning,'','')[3:-1])
                        PARAMS['warnmx'] -= 1
                        if PARAMS['warnmx'] == 0: print('skipping more warnings ...')
    def dispersion(self,steps=10,closed=True):
        """ track the dispersion function """
        traj = []
        # column vector with MDIM rows, 1 column
        v_0 = NP.array([0.,0.,0.,0.,0.,1.,0.,0.,0.,0.])
        if closed == True:
            m_cell = self.acc_node
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
                mx = i_element.matrix
                v_0 = NP.dot(mx,v_0)
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
            viseo = element.viseo
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
    def cs_traj(self,steps=1):
        """ track cosine- & sine-like trajectories """
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
        
        """ injektion parameters """
        tkin = self.injection_energy
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
        # INITIAL @ entrance
        # MDIMxMDIM tracking used here
        s   = 0.
        c_0 = NP.array([x1, x1p, y1, y1p, z1, z1p, tkin,1,0,1])  # C
        s_0 = NP.array([x4, x4p, y4, y4p, z4, z4p, tkin,1,0,1])  # S
        # function names
        c_fun = Functions(('s','cx','cxp','cy','cyp','cz','cdp'))
        s_fun = Functions(('s','sx','sxp','sy','syp','sz','sdp'))
        c_fun.append(s,c_0)
        s_fun.append(s,s_0)

        """ loop through lattice """
        # TODO this function needs serious brain checking again. energy over distance? not carrried correctly! needed?
        for element in iter(self):
            if isinstance(element,ELM.MRK): continue     # TODO skip markers, not a permanent solution
            particle = element.particle
            gamma    = particle.gamma
            tkin     = particle.tkin
            slices   = element.make_slices(anz=steps)
            try:
                for i_element in slices:
                    # if i_element.type == 'QFth' : print(i_element.type,i_element.matrix)
                    # if i_element.type == 'QDth' : print(i_element.type,i_element.matrix)
                    s += i_element.length
                    """ begin map """
                    ## COSine_like
                    c_0 = i_element.map(c_0)   
                    cx  = c_0[XKOO]
                    cxp = c_0[XPKOO]
                    cy  = c_0[YKOO]
                    cyp = c_0[YPKOO]
                    # cz  = -c_0[ZKOO]*360./(beta*lamb)            # sigmaz_i --> dPhi [deg]
                    # cdw = c_0[ZPKOO]*(gamma+1.)/gamma*100.       # dp/p --> dW/W [%]
                    cz  = c_0[ZKOO]*1.e3      # z [mm]
                    cdp = c_0[ZPKOO]*100.     # dp/p [%]
                    c_1 = NP.array([cx, cxp, cy, cyp, cz, cdp, tkin,1,0,1])  # C

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
                    s_1 = NP.array([sx, sxp, sy, syp, sz, sdp, tkin,1,0,1])  # S

                    # c_fun.append(s,(cx,cxp,cy,cyp,cz,cdp))
                    # s_fun.append(s,(sx,sxp,sy,syp,sz,sdp))
                    c_fun.append(s,c_1)
                    s_fun.append(s,s_1)
                    """ end map """
                    if 0: DEBUG_TRACKs(i_element,(cx,cxp,cy,cyp,cz,cdp),(sx,sxp,sy,syp,sz,sdp))
            except (ValueError,ELM.OutOfRadialBoundEx) as ex:
                reason = ex.__class__.__name__
                print('@map in cs_traj: {} at s={:6.2f} [m]'.format(reason,s))
                sys.exit(1)
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
        s = NP.dot(self.acc_node.matrix.T,s)
        s = NP.dot(s,self.acc_node.matrix)
        # dets = LA.det(s)
        # if fabs(dets-1.) > 1.e-12:
            # for i in range(ELM.Matrix._dim):
                # print('[{:4>+.6f}, {:4>+.6f}, {:4>+.6f}, {:4>+.6f}, {:4>+.6f}, {:4>+.6f}]\n'.
                    # format(s[i,0],s[i,1],s[i,2],s[i,3],s[i,4],s[i,5]),end='')
        res = [s[0,1],s[1,0],s[2,3],s[3,2],s[4,5],s[5,4]]
        return(res)
    def concat(self,lattice):
        """Concatenate two Lattice pieces (self+lattice)"""
        for element in iter(lattice):
            self.add_node(element)
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
class TestLattice(unittest.TestCase):
    def test_lattice_add_first_6_nodes(self):
        print('----------------------------------test_lattice_add_first_6_nodes')
        rfgap = ELM.GAP("GAP", 2., -25., 0.044, 800.e6,particle=Proton(66.))
        drift = ELM.D("DR",length=1.)
        quadf = ELM.QF("QUADF",3.,particle=Proton(66.),length=0.5)
        lattice = Lattice(injection_energy=33.)
        lattice.add_node(rfgap)
        lattice.add_node(drift)
        lattice.add_node(quadf)
        lattice.add_node(rfgap)
        lattice.add_node(drift)
        lattice.add_node(quadf)
        for i in range(len(lattice.seq)):
            print(lattice.seq[i].ref_track, lattice.seq[i].particle.tkin,lattice.seq[i].label)
    def test_lattice_generator(self):
        print('----------------------------------test_lattice_generator')
        input_file = "unittests/simuIN.yml"
        lattice = factory(input_file)
        print(lattice.__class__.__name__,F"{len(lattice.seq)} elements")
        print(lattice.__class__.__name__,F"length {len(lattice.seq)} [m]")
        self.assertEqual(lattice.__class__.__name__,"Lattice")
    def test_lattice_iterator(self):
        input_file = "unittests/simuIN.yml"
        print('---------------------------------test_lattice_iterator')
        lattice = factory(input_file)
        # node = None
        for cnt,node in enumerate(iter(lattice)): 
            print(F"left-->right {cnt}\r",end="")
            # print(cnt,'{:38s} {:38s}'.format(repr(node),repr(node.next)))
            pass
        print()
        self.assertEqual(node,lattice.seq[-1],"last node")

        lattice.toggle_iteration()
        for cnt,node in enumerate(iter(lattice)):
            print(F"right-->left {cnt}\r",end="")
            # print(cnt,'{:38s} {:38s}'.format(repr(node),repr(node.prev)))
            pass
        print()
        self.assertEqual(node,lattice.seq[0],"first node")
        self.assertNotEqual(lattice.seq[0],lattice.seq[-1],"first node .ne. last node")
    def test_lattice_concat(self):
        print('---------------------------------test_lattice_concat')
        l = 1.
        p = Proton(50.)
        lattice1 = Lattice()
        for i in range(5):
            lattice1.add_node(ELM.D(F"Drift{i}",length=l,aperture=5.))
        for i in range(5):
            print(lattice1.seq[i].label,lattice1.seq[i].position)
        print()
        lattice2 = Lattice()
        for i in range(5):
            lattice2.add_node(ELM.D(F"Drift{6+i}",length=l,aperture=5.))
        for i in range(5):
            print(lattice2.seq[i].label,lattice1.seq[i].position)
        print()
        lattice1.concat(lattice2)
        for i in range(4,len(lattice1.seq)):
            print(lattice1.seq[i].label,lattice1.seq[i].position)
    def test_first_last_gap(self):
        print('----------------------------------test_first_last_gap')
        input_file = "unittests/simuIN.yml"
        lattice = factory(input_file)
        print(F"first gap {lattice.first_gap.label} {lattice.first_gap.position}")
        lattice.toggle_iteration()
        print(F"last  gap {lattice.first_gap.label} {lattice.first_gap.position}")
    def test_wille_lattice_with_concat(self):
        from matplotlib.pyplot import plot,show,legend
        def make_wille():
            def make_wille_lattice():
                """ Wille's test lattice """
                PARAMS['injection_energy'] = 50.
                kqf   = -1.2
                kqd   = -kqf
                lqf   = 0.20
                lqd   = 0.40
                rhob  = 3.8197
                lb    = 1.50
                phib  = 11.25   #]deg]
                ld    = 0.55
                p     = Proton(PARAMS['injection_energy'])
                gradf = kqf*p.brho
                gradd = kqd*p.brho

                lattice = Lattice()
                anz = 0
                while anz < 1:
                    anz += 1
                    """ nodes """
                    # wedge= ELM.Wedge(phib,rhob,t3d_wedge=False)
                    w1   = ELM.Wedge(phib,rhob,t3d_wedge=False)
                    w2   = ELM.Wedge(phib,rhob,t3d_wedge=False)
                    w3   = ELM.Wedge(phib,rhob,t3d_wedge=False)
                    w4   = ELM.Wedge(phib,rhob,t3d_wedge=False)
                    qf1  = ELM.QF('QF1',gradf,particle=p,length=lqf)
                    qf2  = ELM.QF('QF2',gradf,particle=p,length=lqf)
                    qd1  = ELM.QD('QD1',gradd,particle=p,length=lqd)
                    d1   = ELM.D('D1',particle=p,length=ld)
                    d2   = ELM.D('D2',particle=p,length=ld)
                    d3   = ELM.D('D3',particle=p,length=ld)
                    d4   = ELM.D('D4',particle=p,length=ld)
                    #TODO test rectangular dipoles (ELM>RD) not done yet
                    # br1  = ELM.RD('RD1',2.*phib,rhob,wedge,particle=p)  
                    # br2  = ELM.RD('RD2',2.*phib,rhob,wedge,particle=p)
                    sd1  = ELM.SD('SD1',2.*phib,rhob,particle=p)
                    sd2  = ELM.SD('SD2',2.*phib,rhob,particle=p)

                    # Mz = qf1*d1*w1*sd1*w2*d2*qd1*d3*w3*sd2*w4*d4*qf2      # complete lattice as single Node

                    lattice.add_node(qf1) 
                    lattice.add_node(d1)
                    # lattice.add_node(br1)
                    lattice.add_node(w1)
                    lattice.add_node(sd1)
                    lattice.add_node(w2)
                    lattice.add_node(d2)
                    lattice.add_node(qd1)
                    lattice.add_node(d3)
                    # lattice.add_node(br2)
                    lattice.add_node(w3)
                    lattice.add_node(sd2)
                    lattice.add_node(w4)
                    lattice.add_node(d4)
                    lattice.add_node(qf2)
                return lattice
            lattice  = make_wille_lattice()
            latticeA = make_wille_lattice()
            lattice.concat(latticeA)
            return lattice
        print('----------------------------------test_wille_lattice_with_concat()')
        print("K.Wille's Beispiel auf pp.113 Formel (3.200)")
        PARAMS['emitx_i'] = PARAMS['emity_i'] = 1.e-6
        lattice = make_wille()
        # cell boundaries
        full_cell   = lattice.cell(closed=True)
        betax,a,g,e = PARAMS['twiss_x_i']()
        betay,a,g,e = PARAMS['twiss_y_i']()
        PARAMS['twiss_z_i'] = Twiss(1.,0.,1.)
        lattice.symplecticity()
        # twiss functions
        twiss_fun = lattice.twiss_funcs(steps=5)
        # cl,sl = lattice.cs_traj(steps=100)sK
        disp = lattice.dispersion(steps=5,closed=True)
        # plots
        s  = [twiss_fun(i,'s')  for i in range(twiss_fun.nbpoints)]
        xs = [twiss_fun(i,'bx') for i in range(twiss_fun.nbpoints)]
        ys = [twiss_fun(i,'by') for i in range(twiss_fun.nbpoints)]
        sdisp = [x[0] for x in disp]             # abszisse s
        disp  = [x[1] for x in disp]             # dispersion(s)
        #-------------------- lattice viseo
        lat_plot, ape_plot = lattice.lattice_plot_functions()
        vis_abszisse = [lat_plot(i,'s')              for i in range(lat_plot.nbpoints)]
        vis_ordinate = [lat_plot(i,'viseo')          for i in range(lat_plot.nbpoints)]
        vzero        = [0.                           for i in range(lat_plot.nbpoints)] # zero line

        plot(s,xs,label='betax')
        plot(s,ys,label='betay')
        plot(sdisp,disp,label='disp')
        plot(vis_abszisse,vis_ordinate,label='',color='black')
        plot(vis_abszisse,vzero,color='black')
        legend(loc='upper left')
        show()
if __name__ == '__main__':
    from lattice_generator import factory
    unittest.main()