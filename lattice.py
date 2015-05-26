#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from math import sqrt,fabs,acos,pi,degrees
import numpy as NP
from numpy import linalg as LA
from copy import copy
from pylab import plot,show,legend
from setup import wille,CONF,SUMMARY,Beam,objprnt,printv
import elements as ELM
import warnings

class Lattice(object):
    def __init__(self):
    ## the Lattice object is a sequence of tuples: (ELM.<element>, from_position, to_position)
        self.seq    = [] 
        self.length = 0.
        self.betax0 = 0.
        self.alfax0 = 0.
        self.gammx0 = 0.
        self.betay0 = 0.
        self.alfay0 = 0.
        self.gammy0 = 0.
        self.full_cell = 0.
    def add_element(self,elment):  ## add element to lattice
        if len(self.seq) == 0:
            s0=0.
        else:
            s0=self.seq[-1][-1]
        l=elment.length
        self.length += l
        elm_with_position = (elment,s0,s0+l)
        self.seq.append(elm_with_position)
    def out(self):                 ## simple log lattice layout to stdout (could be better!)
        mcell=ELM.I(label='<=Lattice')     ##  chain matrices 
        for ipos in self.seq:
            element,s0,s1 = ipos
            printv(1,'{:s}\tlength={:.3f}\tfrom-to: {:.3f} - {:.3f}'.
                  format(element.label,element.length,s0,s1))
            mcell = element * mcell   ## Achtung: Reihenfolge im Produkt ist wichtig! Umgekehrt == Blödsinn
        mcell.out()
    def energy_trim(self):         ## trim lattice matrices for changing beam energy
        def min(a,b):
            r = b
            if a < b:
                r = a
            return r
        def max(a,b):
            r = b
            if a > b:
                r = a
            return r
        cav_counter = 0
        qf_counter  = 0
        qd_counter  = 0
        ttfm = +1.e+50
        ttfx = +1.e-50
        seq_trimmed = []
        printv(1,'Beam @ entrance:\n'+Beam.soll.out(tee=False))
        tk_i = Beam.soll.tkin
        for item in self.seq:
            element,s0,s1 = item
            updated = element.update()
            # print(updated.label,'\t',Beam.soll.tkin)
            seq_trimmed.append((updated,s0,s1))
            if isinstance(updated,ELM.QF) and (not isinstance(updated,ELM.QD)):
                qf_counter += 1
            if isinstance(updated,ELM.QD):
                qd_counter += 1
            if isinstance(updated,ELM.RFG) or isinstance(updated,ELM.RFC):
                cav_counter += 1
                ttfm = min(updated.tr,ttfm)
                ttfx = max(updated.tr,ttfx)
        printv(1,'Beam @ exit:\n'+Beam.soll.out(tee=False))
        tk_f = Beam.soll.tkin
        SUMMARY['nboff F-quads']        = qf_counter
        SUMMARY['nboff D-quads']        = qd_counter
        SUMMARY['nboff cavities']       = cav_counter
        SUMMARY['(ttf_min,ttf_max)']    = (ttfm,ttfx)
        SUMMARY['(energy_i,energy_f) [MeV]']  = (tk_i,tk_f)
        self.seq = seq_trimmed      ## replace myself
    def cell(self,closed=True):    ## full cell inspection
        if self.full_cell == 0.0:
            mcell=ELM.I(label=' <= Entrance')     ##  chain matrices for full cell
            for count, ipos in enumerate(self.seq):
                element,s0,s1 = ipos
                mcell = element * mcell   ## Achtung: Reihenfolge im Produkt ist wichtig! Umgekehrt == Blödsinn

            # Stabilität ?
            unstable=False
            stab = fabs(mcell.tracex())
            # if verbose:
            printv(0,'stability X? ',stab)
            if stab >= 2.0:
                # if verbose:
                printv(0,'unstable Lattice in x-plane\n')
                unstable=True
            else:
                cos_mux = 0.5 * stab
                mux = degrees(acos(cos_mux))

            stab = fabs(mcell.tracey())
            # if verbose:
            printv(0,'stability Y? ',stab)
            if stab >= 2.0:
                # if verbose:
                printv(0,'unstable Lattice in y-plane\n')
                unstable=True
            else:
                cos_muy = 0.5 * stab
                muy = degrees(acos(cos_muy))
            if not unstable:
                # if verbose:
                printv(0,'\nphase_advance: X[deg]={:3f} Y[deg]={:.3f}\n'.format(mux,muy))

            self.full_cell = mcell    # the full cell
            # if verbose:
            printv(1,'Lattice.cell: Zellenmatrix (i)->(f)')
            self.full_cell.out()
            det = LA.det(self.full_cell.matrix)
            # if verbose:
            printv(1,'det|full-cell|={:.5f}\n'.format(det))
            # Determinate M-I == 0 ?
            beta_matrix = mcell.BetaMatrix()
            for i in range(5):
                beta_matrix[i,i] = beta_matrix[i,i]-1.0
            det = LA.det(beta_matrix)
            # if verbose:
            printv(1,'det|Mbeta - I|={:.5f}\n'.format(det))
            # symplectic?
            s = self.symplecticity()
            # if verbose:
            printv(1,'symplectic (+1,-1,+1,-1,+1,-1)?')
            printv(1,'[{:4>+.2f}, {:4>+.2f}, {:4>+.2f}, {:4>+.2f}, {:4>+.2f}, {:4>+.2f}]\n'.
                format(s[0],s[1],s[2],s[3],s[4],s[5]))
            
            # Startwerte für twiss-functions aus Eigenwert- und Eigenvektor
            # beta_matrix = mcell.BetaMatrix()  
            # eigen, vectors = LA.eig(beta_matrix)
            # print('Eigenwerte\n',eigen)
            # print('Eigenvektoren\n',vectors)
            # show them !
            # for i in (0,3):
                # print('Mit numpy.linalg berechneter Eigenwert: \n',eigen[i].real)
                # bx=vectors[0][i].real; ax=vectors[1][i].real; gx=vectors[2][i].real
                # by=vectors[3][i].real; ay=vectors[4][i].real; gy=vectors[5,i].real
                # print('...und sein Eigenvektor dazu: \n {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}'.
                # format(bx,ax,gx,by,ay,gy))
            # probe: M*beta = beta for full cell ?
            # print(vectors[:,i].T)
            # probe=beta_matrix.dot(vectors[:,i].T)  
            # print(probe)  
            # Linearkombination beider Eigenvektoren zum Eigenwert 1 
            # vector = vectors[:,0]+vectors[:,3]
            # bax = vector[0].real
            # Vorzeichenkorrektur
            # sign_betax = -1.0 if bax <= 0.0 else +1.0
            # bay = vector[3].real
            # sign_betay = -1.0 if bay <= 0.0 else +1.0
            # neue Linearkombination damit twiss-functions positiv sind
            # vector = sign_betax*vectors[:,0]+sign_betay*vectors[:,3]
            # print('Summe beider Eigenvektoren zum reellen Eigenwert\n',vector)

            # das sollten die twiss-functions am Ein- and Ausgang sein...
            # bax = vector[0].real
            # alx = vector[1].real
            # gmx = vector[2].real
            # bay = vector[3].real
            # aly = vector[4].real
            # gmy = vector[5].real

            if closed:
                if not unstable:
                    # Startwerte für twiss-functions aus Formeln von K.Wille (Teubner Studienbücher)
                    cell_matrix = self.full_cell.matrix
                    m11 =cell_matrix[0,0];  m12 =cell_matrix[0,1]
                    m21 =cell_matrix[1,0];  m22 =cell_matrix[1,1]
                    n11 =cell_matrix[2,2];  n12 =cell_matrix[2,3]
                    n21 =cell_matrix[3,2];  n22 =cell_matrix[3,3]
                    bax=fabs(2.0-m11*m11-2.*m12*m21-m22*m22)
                    bax=sqrt(bax)
                    bax=2.0*m12/bax
                    if(bax < 0.):
                        bax= -bax
                    alx=(m11-m22)/(2.*m12)*bax
                    gmx=(1.+alx*alx)/bax
                    bay=fabs(2.0-n11*n11-2.*n12*n21-n22*n22)
                    bay=sqrt(bay)
                    bay=2.0*n12/bay
                    if(bay < 0.):
                        bay = -bay
                    aly=(n11-n22)/(2.*n12)*bay
                    gmy=(1.+aly*aly)/bay
                    
                    # Probe: twiss-functions durch ganze Zelle (nur sinnvoll fuer period. Struktur!)
                    v_beta=NP.array([[bax],[alx],[gmx],[bay],[aly],[gmy]])
                    m_cell=self.full_cell.BetaMatrix()
                    v_beta_end = m_cell.dot(v_beta)
                    # if verbose:
                    printv(0,'Probe: {Twiss_Ende} == {Zellenmatrix}x{Twiss_Anfang}?')
                    printv(0,'Anfang: ',v_beta.T)
                    printv(0,'Ende  : ',v_beta_end.T,'\n')
                    CONF['sigx_i'] = sqrt(bax*CONF['emitx_i'])
                    CONF['sigy_i'] = sqrt(bay*CONF['emity_i'])
                else:
                    raise RuntimeError('stop: unstable lattice')
            else:
                # Startwerte fuer transfer line (keine periodischen Randbedingungen!)
                xi=CONF['sigx_i']
                yi=CONF['sigy_i']
                alx=aly=0.
                emix = CONF['emitx_i']
                emiy = CONF['emity_i']
                bax=xi**2/emix
                bay=yi**2/emiy
                gmx=(1.+alx*alx)/bax
                gmy=(1.+aly*aly)/bay
                
            self.betax0 = bax
            self.alfax0 = alx
            self.gammx0 = gmx
            self.betay0 = bay
            self.alfay0 = aly
            self.gammy0 = gmy
            
        return (self.full_cell,self.betax0,self.betay0)
    def report(self):              ## lattice layout report (may not work!)
        raise RuntimeWarning('Lattice.report() not ready')
        reprt = ''
        header = ''
        row = ''
        for count, ipos in enumerate(reversed(self.seq)):
            elm,si,sf = ipos
            name= elm.label
            len =elm.length
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
    def reverse(self):             ## return a reversed Lattice (not used! probably bogus!)
        raise RuntimeWarning('Lattice.reverse() not ready')
        res=Lattice()
        seq=copy(self.seq)
        seq.reverse()
        for ipos in seq:
            elm,s,s=ipos
            res.add_element(elm)
        return res
    def append(self,piece):        ## concatenate two Lattice pieces
        seq=copy(piece.seq)  
        for ipos in seq:
            elm,s0,s1=ipos
            self.add_element(elm)
    def functions(self,steps=10):  ## functions of s
        beta_fun=[]
        # ms=ELM.I()
        bx = self.betax0
        ax = self.alfax0
        gx = self.gammx0
        by = self.betay0
        ay = self.alfay0
        gy = self.gammy0
        v_beta0=NP.array([[bx],[ax],[gx],[by],[ay],[gy]])
        v_beta = v_beta0
        # print(v_beta0)
        s=0.0
        for ipos in self.seq:
            element,s0,s1 = ipos
            for count,i_element in enumerate(element.step_through(steps)):
                # ms = i_element*ms
                # print('i_element.label={} viseo={} ms.length={}'.format(i_element.label,i_element.viseo,ms.length))
                # m_beta = ms.BetaMatrix()
                # v_beta = m_beta.dot(v_beta0)
                m_beta = i_element.BetaMatrix()
                v_beta = m_beta.dot(v_beta)
                s += i_element.length
                betax = v_beta[0,0]
                betay = v_beta[3,0]
                viseo = i_element.viseo
                beta_fun.append((s,betax,betay,viseo))
        (c_like,s_like) = self.cs_traj(steps)
        return (beta_fun,c_like,s_like)
    def dispersion(self,steps=10,closed=True):  ## dispersion (not used! probably bogus!)
        warnings.warn(UserWarning('Lattice.dispersion() not ready'))
        traj=[]
        v_0=NP.array([[0.],[0.],[0.],[0.],[0.],[1.]])
        if closed == True:
            m_cell = self.full_cell
            m11=m_cell.matrix[0,0]
            m15=m_cell.matrix[0,5]
            d0 = m15/(1.-m11)     # from H.Wiedemann (6.79) pp.206
            v_0=NP.array([[d0],[0.],[0.],[0.],[0.],[1.]])
        s=0.0
        for ipos in self.seq:
            element,s0,s1 = ipos
            for i_element in element.step_through(steps):
                m_beta = i_element.matrix
                v_0 = m_beta.dot(v_0)
                s += i_element.length
                d  = v_0[0,0]
                dp = v_0[1,0]
                viseo = i_element.viseo
                traj.append((s,d,dp))
        return traj
    def cs_traj(self,steps=10):    ## cosine & sine trajectories 
        lamb = CONF['wellenlänge']
        c_like =[]
        s_like =[]
        x1  = sqrt(CONF['emitx_i']*self.betax0) # x-plane: principal-1 (cos like) with alpha=0
        x1p = 0.
        x2  = 0.                                # y-plane: principal-2 (sin like)
        x2p = sqrt(CONF['emitx_i']/self.betax0)
        y1  = sqrt(CONF['emity_i']*self.betay0)
        y1p = 0.
        y2  = 0.
        y2p = sqrt(CONF['emity_i']/self.betay0)
        dz  = CONF['dZ']                        # eingabe dZ
        dp  = CONF['dP/P']                      # eingabe dP/P0
        c_0=NP.array([[x1],[x1p],[y1],[y1p],[dz],[0.]])     # cos-like traj.
        s_0=NP.array([[x2],[x2p],[y2],[y2p],[0.],[dp]])     # sin-like traj.
        s=0.0
        for ipos in self.seq:
            element,s0,s1 = ipos
            beam = element.beam
            # objprnt(beam,text=element.label)
            for i_element in element.step_through(steps):
                element_matrix = i_element.matrix
                c_0 = element_matrix.dot(c_0)
                s_0 = element_matrix.dot(s_0)
                s += i_element.length
                # cos_like
                cx  = c_0[0,0]
                cxp = c_0[1,0]
                cy  = c_0[2,0]
                cyp = c_0[3,0]
                cz  = -c_0[4,0]*360./(beam.beta*lamb)           # conversion zu dPhi [deg]
                cdw = c_0[5,0]*(beam.gamma+1.)/beam.gamma*100.  # conversion zu dW/W [%]
                # sin_like
                sx  = s_0[0,0]
                sxp = s_0[1,0]
                sy  = s_0[2,0]
                syp = s_0[3,0]
                sz  = -s_0[4,0]*360./(beam.beta*lamb)
                sdw = s_0[5,0]*(beam.gamma+1.)/beam.gamma*100.
                c_like.append((cx,cxp,cy,cyp,cz,cdw))
                s_like.append((sx,sxp,sy,syp,sz,sdw))
        return (c_like,s_like)
    def symplecticity(self):       ## test symplecticity
        s=NP.array([[0., 1., 0.,0., 0.,0.],
                    [-1.,0., 0.,0., 0.,0.],
                    [ 0.,0., 0.,1., 0.,0.],
                    [ 0.,0.,-1.,0., 0.,0.],
                    [ 0.,0., 0.,0., 0.,1.],
                    [ 0.,0., 0.,0.,-1.,0.]])
        s=NP.dot(self.full_cell.matrix.T,s)
        s=NP.dot(s,self.full_cell.matrix)
        # dets=LA.det(s)
        # if fabs(dets-1.) > 1.e-12:
            # for i in range(ELM.Matrix._dim):
                # print('[{:4>+.6f}, {:4>+.6f}, {:4>+.6f}, {:4>+.6f}, {:4>+.6f}, {:4>+.6f}]\n'.
                    # format(s[i,0],s[i,1],s[i,2],s[i,3],s[i,4],s[i,5]),end='')
        res=[s[0,1],s[1,0],s[2,3],s[3,2],s[4,5],s[5,4]]
        return(res)
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
def make_lattice():  # a test lattice
     print("K.Wille's Beispiel auf pp. 112-113")
     kqf=  wille()['k_quad_f']
     lqf=  wille()['length_quad_f']
     kqd=  wille()['k_quad_d']
     lqd=  wille()['length_quad_d']
     rhob= wille()['beding_radius'] 
     lb=   wille()['dipole_length']
     ld=   wille()['drift_length']
     ## elements
     mqf=ELM.QF(kqf,lqf,'QF')
     mqd=ELM.QD(kqd,lqd,'QD')
     mb=ELM.SD(rhob,lb,'B')
     mb1=ELM.SD(rhob,lb*0.5,'B1')  ## 1/2 sector dip.
     mw=ELM.WD(mb)
     mw1=ELM.WD(mb1)
     mbr=ELM.RD(rhob,lb)
     md=ELM.D(ld)    
     ## lattice
     lattice=Lattice()
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
     # lattice.out()
     top=Lattice()
     top.append(lattice)
     # top.append(top)
     # top.append(top)
     # top.append(top)
     # top.append(top)
     # top.append(top)
     # top.append(top)
     return top
def test0():
    lat = make_lattice()
    lat.out()
    print('--------------- EOF test0 --------------------')
def test1():
    lattice=make_lattice()
    mcell,betax,betay=lattice.cell()
    beta_matrix = mcell.BetaMatrix()
    
    eigen, vectors = LA.eig(beta_matrix)
    print('eigen\n',eigen)
    print('vectors\n',vectors)
    print('Mit numpy.linalg berechneter Eigenwert: \n',eigen[0].real)
    bx=vectors[0][0].real; ax=vectors[1][0].real; gx=vectors[2][0].real
    print('...und sein Eigenvektor dazu: \n {:.6f} {:.6f} {:.6f}'.format(bx,ax,gx))
    probe=beta_matrix.dot(vectors[:,0])
    print('Probe: \n {:.6f} {:.6f} {:.6f}'.format(probe[0].real,probe[1].real,probe[2].real))
    # print(
# '''Ergebnis scheint trivial zu sein!
# Eigenvvektor ist auf 1 normiert.
# Jedes Vielfache eines Eigenvektors ist auch Eigenvektor!
# Wo ist der Wert von betax geblieben?
# Antwort: Da die Emittanz noch nicht bekannt ist kann man den
# Eigenvektor auf 1 normieren und erhält dasselbe wie mit Wille's Formeln.''')
    print('--------------- EOF test1 --------------------')
def test2():
    lattice=make_lattice()
    ## cell boundaries
    mcell,betax,betay=lattice.cell(closed=True)
    print('BETAx[0] {:.3f} BETAy[0] {:.3f}'.format(betax,betay))
    lattice.symplecticity()
    ## lattice function as f(s)
    beta_fun,cl,sl = lattice.functions(steps=100)   
    disp = lattice.dispersion(steps=100,closed=True)
    ## plots
    s  = [x[0] for x in beta_fun]    # s
    xs = [x[1] for x in beta_fun]    # betax
    ys = [x[2] for x in beta_fun]    # betay
    ds = [x[1] for x in disp]        # dispersion
    vs = [x[3]-1. for x in beta_fun] # viseo
    zero=[-1. for x in beta_fun]     # viseo
##    for i in range(len(s)):
##        print('s, betax(s) betay(s)',s[i],xs[i],ys[i])
    plot(s,xs,label='bx/bx0')
    plot(s,ys,label='by/by0')
    plot(s,ds,label='dp/p')
    plot(s,vs,label='element',color='black')
    plot(s,zero,color='black')
    legend(loc='upper left')
    show(block=True)
    print('--------------- EOF test2 --------------------')
if __name__ == '__main__':
    test0()
    test1()
    test2()
