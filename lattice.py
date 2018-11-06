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
from math import sqrt,fabs,acos,degrees
from numpy import linalg as LA
import numpy as NP
from copy import copy
import warnings

from setutil import XKOO, XPKOO, YKOO, YPKOO, ZKOO, ZPKOO, EKOO, DEKOO, SKOO, LKOO
from setutil import wille,PARAMS,FLAGS,SUMMARY,printv,DEBUG,sigmas, objprnt, Ktw, Ktp, Twiss
import elements as ELM
import TTFG as TTF
from sigma import Sigma

## DEBUGING
def DEBUG_ON(*args):
    DEBUG(*args)
def DEBUG_OFF(*args):
    pass
DEBUG_MODULE = DEBUG_OFF

#todo: cos-like, sin-like traj in z ???      
## Lattice
class Lattice(object):
    """
    The Lattice object is a list of elements: ELM.<element> in self.seq
    """
    def __init__(self):
        self.seq    = []
        self.length = 0.
        self.accel  = 0.
    
    def __getitem__(self,n):
        item = self.seq[n]
        return item

    def add_element(self,element):
        """ add element to lattice """
        if len(self.seq) == 0:
            s0 = 0.
        else:
            s0 = self.seq[-1].position[2]
        l = element.length
        si = s0
        sf = si+l
        sm = (sf+si)/2.
        position = (si,sm,sf)
        element.position = position 
        self.length = sf
        # DEBUG('add_element: ',' [si,sm,sf,]=[{1}] {0}'.format(repr(element),''.join('{:5.3f},'.format(el) for el in element.position)))
        self.seq.append(element)

    def string(self):
        """ log lattice layout to string (could be even better?) """
        mcell = ELM.I(label='')   #  chain matrices
        for element in self.seq:
            DEBUG_MODULE('{:10s}({:d})\tlength={:.3f}\tfrom-to: {:.3f} - {:.3f}'.format(element.label,id(element),element.length,element.position[0],element.position[2]))
            mcell = element * mcell   # Achtung: Reihenfolge im Produkt ist wichtig! Umgekehrt == Blödsinn
        mcell.section = '<= full lattice map'
        return mcell.string()

    def stats(self,soll_track):
        """ gather lattice statistics """
        cav_counter = 0
        q_counter   = 0
        ttfm = +1.e+50
        ttfx = +1.e-50
        tk_i = soll_track.getpoints()[0]()[6]
        tk_f = soll_track.getpoints()[-1]()[6]
        for element in self.seq:
            if isinstance(element,(ELM.QF,ELM.QD)):
                q_counter += 1
            if isinstance(element,(ELM.RFG,ELM.RFC)):
                cav_counter += 1
                ttfm = min(element.tr,ttfm)
                ttfx = max(element.tr,ttfx)
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
        for count,element in enumerate(self.seq):
            # Achtung: Reihenfolge im Produkt ist wichtig! Umgekehrt == Blödsinn
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
        emitx = PARAMS['emitx_i']
        emity = PARAMS['emity_i']

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
                print('unstable lattice - STOP')
                sys.exit(1)
        else:
            ## transversale twiss parameter fuer transfer lines
            # alfa, beta und emittance definieren den beam @ entrance, NOTE: transfer lattices need not to be stable!
            bax,alx,gmx,epsx = PARAMS['twiss_x_i']()
            bay,aly,gmy,epsy = PARAMS['twiss_y_i']()
        printv(0,'using @ entrance: [beta,  alfa,  gamma]-X    [beta,   alfa,   gamma]-Y')
        printv(0,'                  [{:.3f}, {:.3f}, {:.3f}]-X    [{:.3f},  {:.3f},  {:.3f}]-Y'.format(bax,alx,gmx,bay,aly,gmy))
        return(self)

    def report(self):
        """ report lattice layout (may not work!) """
        raise RuntimeWarning('Lattice.report() not ready')
        reprt = ''
        header = ''
        row = ''
        for count,element in enumerate(reversed(self.seq)):
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

    def reverse(self):
        raise RuntimeWarning('Lattice.reverse() not implemented and not used! (probably bogus!)')
        res = Lattice()
        seq = copy(self.seq)
        seq.reverse()
        for elm in seq:
            res.add_element(elm)
        return res

    def concat(self,lattice_piece):
        """Concatenate two Lattice pieces"""
        for element in lattice_piece.seq:
            element = copy(element) if element in self.seq else element
            self.add_element(element)

    def twiss_envelopes(self,steps=1):
        """
        Calulate envelopes from initial twiss-vector with beta-matrices
        """
        bx,ax,gx,epsx = PARAMS['twiss_x_i']()
        by,ay,gy,epsy = PARAMS['twiss_y_i']()
        bz,az,gz,epsz = PARAMS['twiss_z_i']()
        twv0 = NP.array([bx,ax,gx,by,ay,gy,bz,az,gz])   # initial
        beta_fun = [list(twv0)+[0]]
        B_matrix = NP.eye(9)                            # cumulated beta-matrix
        for node in iter(self.seq):
            slices = node.make_slices(anz = steps)
            means = []
            s,sm,sf = node.position
            for slice in slices:
                cx  = slice.matrix[Ktp.x ,Ktp.x];    sx  = slice.matrix[Ktp.x ,Ktp.xp]
                cxp = slice.matrix[Ktp.xp,Ktp.x];    sxp = slice.matrix[Ktp.xp,Ktp.xp]
                cy  = slice.matrix[Ktp.y ,Ktp.y];    sy  = slice.matrix[Ktp.y ,Ktp.yp]
                cyp = slice.matrix[Ktp.yp,Ktp.y];    syp = slice.matrix[Ktp.yp,Ktp.yp]
                cz  = slice.matrix[Ktp.z ,Ktp.z];    sz  = slice.matrix[Ktp.z ,Ktp.zp]
                czp = slice.matrix[Ktp.zp,Ktp.z];    szp = slice.matrix[Ktp.zp,Ktp.zp]
    
                beta_matrix = NP.eye(9)
                beta_matrix[0,0] = cx**2;    beta_matrix[0,1] = -2*cx*sx;      beta_matrix[0,2] = sx**2
                beta_matrix[1,0] = -cx*cxp;  beta_matrix[1,1] = sx*cxp+sxp*cx; beta_matrix[1,2] = -sx*sxp
                beta_matrix[2,0] = cxp**2;   beta_matrix[2,1] = -2*cxp*sxp;    beta_matrix[2,2] = sxp**2
    
                beta_matrix[3,3] = cy**2;    beta_matrix[3,4] = -2*cy*sy;      beta_matrix[3,5] = sy**2
                beta_matrix[4,3] = -cy*cyp;  beta_matrix[4,4] = sy*cyp+syp*cy; beta_matrix[4,5] = -sy*syp
                beta_matrix[5,3] = cyp**2;   beta_matrix[5,4] = -2*cyp*syp;    beta_matrix[5,5] = syp**2
                
                B_matrix = NP.dot(beta_matrix,B_matrix)
                twv = NP.dot(B_matrix,twv0)             # track twiss-vector
                s += slice.length
                beta_fun.append(list(twv)+[s])
                bx = twv[Ktw.bx]; ax = twv[Ktw.ax]; gx = twv[Ktw.gx]
                by = twv[Ktw.by]; ay = twv[Ktw.ay]; gy = twv[Ktw.gy]
                sigxy = (*sigmas(ax,bx,epsx),*sigmas(ay,by,epsy))
                means.append(sigxy)
            means = NP.array(means)
            means = NP.mean(means,axis=0)
            node['sigxy'] = means                       # each node has its average sigmas
            # aperture check
            if FLAGS['useaper']:
                n_sigma = PARAMS['n_sigma']
                if node.aperture != None:
                    aperture = node.aperture
                    sigx, sigxp, sigy, sigyp = node['sigxy']
                    if(aperture < n_sigma*sigx or aperture < n_sigma*sigy):
                        warnings.showwarning(
                            '{} sigma aperture hit @ s={:.1f} [m]'.format(n_sigma,sm),
                            UserWarning,'lattice.py',
                            'twiss_functions()')
        return beta_fun
        
    def sigmas(self,steps = 10):
        """ dispatch to different envelope functions """
        def envelopes(function, steps = 10):
            """ calc. envelopes using function """
            beta_fun  = function(steps = steps)
            b,a,g,epsx = PARAMS['twiss_x_i']()
            b,a,g,epsy = PARAMS['twiss_y_i']()
            # for x in beta_fun:
            #     print(x[Ktw.bx])
            sigma_fun = [(x[Ktw.s],sqrt(x[Ktw.bx]*epsx),sqrt(x[Ktw.by]*epsy)) for x in beta_fun]
            return sigma_fun

        if FLAGS['sigma']:
            mess = 'CALCULATE SIGMA ENVELOPES'
            function = self.sigma_envelopes # use sigma-matrix
        elif not FLAGS['sigma']:
            mess = 'CALCULATE TWISS ENVELOPES'
            # function = self.twiss_functions # use beta-matrix
            function = self.twiss_envelopes # use beta-matrix

        if not FLAGS['KVout']: 
            print(mess)
            sigma_fun = envelopes(function, steps = steps)
        return sigma_fun

    def sigma_envelopes(self, steps = 1):
        """ track the sigma-matrix through the lattice and extract twiss functions """
        # initials
        bx,ax,gx,epsx = PARAMS['twiss_x_i']()
        by,ay,gy,epsy = PARAMS['twiss_y_i']()
        bz,az,gz,epsz = PARAMS['twiss_z_i']()
        #todo: check initial values
        twv0     = NP.array([bx,ax,gx,by,ay,gy,bz,az,gz])  # twiss vector IN lattice
        sg0      = Sigma(twv0,epsx,epsy,epsz)              # sigma object IN lattice
        # sigma envelopes as function of distance s
        sigma_fun = []
        for node in self.seq: # loop nodes
            # sigma-matrices for a single node
            sigmas = node.sigma_beam(steps = steps, sg = sg0) 
            # prep plot list of ftn's
            for sg,s in sigmas:
                v = sg.twiss()      # twiss from Sigma object
                flist = v.tolist()
                flist.append(s)
                sigma_fun.append(flist)
            sg0 = sg          # loop back nodes

            # aperture check
            if FLAGS['useaper']:
                n_sigma = PARAMS['n_sigma']
                if node.aperture != None:
                    aperture = node.aperture
                    sigx, sigxp, sigy, sigyp = node['sigxy']
                    si,sm,sf                 = node.position
                    if(aperture < n_sigma*sigx or aperture < n_sigma*sigy):
                        warnings.showwarning(
                            '{} sigma aperture hit @ s={:.1f} [m]'.format(n_sigma,sm),
                            UserWarning,'lattice.py',
                            'sigma_functions()')
        return sigma_fun

    def dispersion(self,steps=10,closed=True):
        """ track the dispersion function """
        traj = []
        v_0 = NP.array([0.,0.,0.,0.,0.,1.,0.,0.,0.,0.])    # column vector with MDIM rows, 1 column
        if closed == True:
            m_cell = self.accel
            m11 = m_cell.matrix[0,0]
            m15 = m_cell.matrix[0,5]
            d0  =  m15/(1.-m11)     # from H.Wiedemann (6.79) pp.206
            v_0[0] = d0
        s = 0.0
        traj = [(s,v_0[0],v_0[1])]
        for element in self.seq:
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
        fun = []   # is list((s = Abzisse,f = Ordinate))
        ape = []
        for element in self.seq:
            # DEBUG((element.__class__,element['viseo'],element.position))
            pos   = element.position
            # element plot
            viseo = element['viseo']
            si, sm, sf = pos
            fun.append((si,0))
            fun.append((si,viseo))
            fun.append((sf,viseo))
            fun.append((sf,0))

            # aperture plot
            aperture = None
            if hasattr(element, 'aperture'): 
                aperture = element.aperture
            if element.__class__.__name__ == 'D': continue
            if aperture == None: continue
            ape.append((sm,aperture))
        return fun,ape

    def cs_traj(self,steps=10):
        """ track cos- & sin-trajectories """
        def SollTest_ON(arg):  # set all 0. to simulate Sollteilchen
            return 0.
        def SollTest_OFF(arg):
            return arg
        soll_test   = SollTest_OFF

        print('CALCULATE C+S TRAJECTORIES')
        gamma       = PARAMS['sollteilchen'].gamma
        beta        = PARAMS['sollteilchen'].beta
        tkin        = PARAMS['sollteilchen'].tkin
        lamb        = PARAMS['wellenlänge']
        
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
        z1,z1p = soll_test((sigmaz_i, 0.))
        z4,z4p = soll_test((0., Dp2p_i))
        # MDIMxMDIM tracking used here
        s      = 0.
        c_like = []; s_like = []
        c_0 = NP.array([x1, x1p, y1, y1p, z1, z1p, tkin,1,0,1])  # C
        s_0 = NP.array([x4, x4p, y4, y4p, z4, z4p, tkin,1,0,1])  # S
        for element in self.seq:
            particle = element.particle
            gamma = particle.gamma
            # objprnt(particle,text='cs_traj: '+element.label)         # DEBUG
            slices = element.make_slices(anz=steps)
            for i_element in slices:
                s += i_element.length
                ## COS_like
                # DEBUG_MODULE('cs_traj: calls {}.map() for C'.format(i_element))
                c_0 = i_element.map(c_0)   # map!!!
                cx  = c_0[XKOO]
                cxp = c_0[XPKOO]
                cy  = c_0[YKOO]
                cyp = c_0[YPKOO]
                cz  = -c_0[ZKOO]*360./(beta*lamb)            # conversion sigmaz_i --> dPhi [deg]
                cdw = c_0[ZPKOO]*(gamma+1.)/gamma*100.       # conversion dp/p --> dW/W [%]
                ## SIN_like
                # DEBUG_MODULE('cs_traj: calls {}.map() for S'.format(i_element))
                s_0 = i_element.map(s_0)   # map!!!
                sx  = s_0[XKOO]
                sxp = s_0[XPKOO]
                sy  = s_0[YKOO]
                syp = s_0[YPKOO]
                sz  = -s_0[ZKOO]*360./(beta*lamb)
                sdw = s_0[ZPKOO]*(gamma+1.)/gamma*100.
                c_like.append((s,cx,cxp,cy,cyp,cz,cdw))
                s_like.append((s,sx,sxp,sy,syp,sz,sdw))
        return (c_like,s_like)

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
    
    @property
    def first_gap(self):
        node = None
        for elm in self.seq:
            if isinstance(elm,(ELM.RFG,ELM.RFC,ELM.GAP)):
                node = elm
                break
        return node

# The commented code is *legacy*. No use to define a new
# subclass and to cast from base class to subclass
# although it worked well!
# 
#     def get_section(self,sec):    *legacy*
#         if not FLAGS['sections']:
#             section = self
#             Section.cast(section)             #the whole lattice is one section
#             setction.set_name('LINAC')
#         else:
#             section = Section(name=sec)
#             for elm in self.seq:
#                 try:
#                     elmsec = elm.get_section()
#                 except AttributeError:
#                     print('WARNING: element {} w/o section attribute. - STOP!'.format(elm.label))
#                     continue
#                 if elmsec == sec:
#                     section.add_element(elm)
#         return section
#     def get_sections(self):
#         sections = []
#         if not FLAGS['sections']:
#             section = self
#             Section.cast(section)             #the whole lattice is one section
#             section.set_name('LINAC')
#             sections.append(section)
#         else:
#             for isec in PARAMS['sections']:
#                 sec = self.get_section(isec)
#                 sections.append(sec)
#         return sections
#        
# class Section(Lattice):
#     """
#     A Lattice with a name
#     """
#     def __init__(self,name='LINAC'):
#         super().__init__()
#         self.name = name
#     def get_name(self):
#         return self.name
#     def set_name(self,name):
#         self.name = name
#     @classmethod
#     def cast(cls,obj):
#         """
#         Convert a BaseClass object into a SubClass object ==> der Trick:
#         ==> cast 'obj' (must be of class Lattice) to object of class Section.
#         """
#         if not isinstance(obj,Lattice):
#             print('ERROR: cast to class Section not possible. -- STOP!')
#             sys.exit(1)
#         obj.__class__ = Section

## Sections
# To add Sections to the lattice I augment the Lattice class with member-functions 
# using the built-in 'setattr(..)'
def get_section(self,sec=None):
    if not FLAGS['sections']:
        section = self       #the whole lattice is one section
        section.name = 'LINAC'
        return section
    else:
        section = Lattice()
        for elm in self.seq:
            try:
                elmsec = elm.section
            except AttributeError as ex:
                print('WARNING: element {} w/o section attribute!'.format(elm.label))
                continue
            if elmsec == sec:
                section.add_element(elm)
    section.name = sec
    return section

def get_sections(self):
    sections = []
    if not FLAGS['sections']:
        sections.append(self.get_section())
    else:
        for isec in PARAMS['sections']:
            sectn = self.get_section(sec=isec)
            sections.append(sectn)
    return sections

#Lattice.get_section  = get_section                 #add method to class Lattice (the wdk way)
#Lattice.get_sections = get_sections                #add method to class Lattice (the wdk way)

setattr(Lattice,get_section.__name__,get_section)   #add method to class Lattice (the python way)
setattr(Lattice,get_sections.__name__,get_sections) #add method to class Lattice (the python way)


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
    mqf1 = ELM.QF(kqf,lqf,'QF1')
    mqf2 = ELM.QF(kqf,lqf,'QF2')
    mqd1 = ELM.QD(kqd,lqd,'QD1')
    md1  = ELM.D(ld)
    md2  = ELM.D(ld)
    md3  = ELM.D(ld)
    md4  = ELM.D(ld)
    mbr1  = ELM.RD(rhob,lb)
    mbr2  = ELM.RD(rhob,lb)
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
#     DEBUG('lattice: ',lattice.string())
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
    s  = [x[Ktw.s]  for x in beta_fun]    # abzisse s
    xs = [x[Ktw.bx] for x in beta_fun]    # betax(s)
    ys = [x[Ktw.by] for x in beta_fun]    # betay(s)
    sd = [x[0] for x in disp]             # abszisse s
    ds = [x[1] for x in disp]             # dispersion(s)
    #-------------------- lattice viseo
    lat_plot, ape_plot = lattice.lattice_plot_functions()
    vsbase = -1.
    vis_abzisse  = [x[0] for x in lat_plot]
    vis_ordinate = [x[1]+vsbase for x in lat_plot]
    vzero        = [vsbase   for x in lat_plot]      # zero line

    plot(s,xs,label='betax')
    plot(s,ys,label='betay')
    plot(sd,ds,label='disp')
    plot(vis_abzisse,vis_ordinate,label='',color='black')
    plot(vis_abzisse,vzero,color='black')
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
if __name__ == '__main__':
    test1()
    test2()
