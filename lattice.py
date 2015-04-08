# -*- coding: utf-8 -*-
import setup as SETUP
import elements as ELM

from math import sqrt, fabs
import numpy as np
from numpy import linalg as LA
from copy import copy
from pylab import plot, show, legend


class Lattice(object):
    def __init__(self):
        self.seq=[] ## sequence of tuples (element, start-position, end-position)
        self.length = 0.
        self.betax0 = 0.
        self.alfax0 = 0.
        self.gammx0 = 0.
        self.betay0 = 0.
        self.alfay0 = 0.
        self.gammy0 = 0.
        self.full_cell = 0.
    #-----------------------
    def add_element(self,elment):
        if len(self.seq) == 0:
            s0=0.
        else:
            s0=self.seq[-1][-1]
        l=elment.length
        self.length += l
        elm_with_position = (elment,s0,s0+l)
        self.seq.append(elm_with_position)
    #-----------------------
    def out(self):
        for ipos in self.seq:
            element,s0,s1 = ipos
            print('{:s} length {:.3f} from {:.3f} to {:.3f}'.
                  format(element.label,element.length,s0,s1))
    #-----------------------
    def cell(self):  ## full cell
        if self.full_cell == 0.0:
            mcell=ELM.I()
            for count, ipos in enumerate(self.seq):
                element,s0,s1 = ipos
                mcell = element * mcell   ## Achtung: Reihenfolge im Produkt ist wichtig! Umgekehrt == Blödsinn

            # Stabilität ?
            stab = fabs(mcell.tracex())
            if stab >= 2.0:
                raise RuntimeError('unstable Lattice in x-plane')
            else:
                print('stability x ',stab)
            stab = fabs(mcell.tracey())
            if stab >= 2.0:
                raise RuntimeError('unstable Lattice in y-plane')
            else:
                print('stability y ',stab)

            # Determinate M-I == 0 ?
            beta_matrix = mcell.BetaMatrix()
            det = LA.det(beta_matrix)
            # print('det {:.5f}\n'.format(det))
            for i in range(5):
                beta_matrix[i,i] = beta_matrix[i,i]-1.0
            det = LA.det(beta_matrix)
            print('det(Mbeta - I) {:.5f}\n'.format(det))
            
            self.full_cell = mcell    # the full cell

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

            ## Startwerte für twiss-functions aus Formeln von K.Wille (Teubner Studienbücher)
            print('Lattice.cell: ganze Zellenmatrix:\n')
            self.full_cell.out()
            cell_matrix = self.full_cell.matrix
            m11 =cell_matrix[0,0];  m12 =cell_matrix[0,1]
            m21 =cell_matrix[1,0];  m22 =cell_matrix[1,1]
            n11 =cell_matrix[2,2];  n12 =cell_matrix[2,3]
            n21 =cell_matrix[3,2];  n22 =cell_matrix[3,3]
            bax=2.0-m11*m11-2.*m12*m21-m22*m22
            bax=sqrt(bax)
            bax=2.0*m12/bax
            if(bax < 0.):
                bax= -bax
            alx=(m11-m22)/(2.*m12)*bax
            gmx=(1.+alx*alx)/bax
            bay=2.0-n11*n11-2.*n12*n21-n22*n22
            bay=sqrt(bay)
            bay=2.0*n12/bay
            if(bay < 0.):
                bay = -bay
            aly=(n11-n22)/(2.*n12)*bay
            gmy=(1.+aly*aly)/bay

            self.betax0 = bax
            self.alfax0 = alx
            self.gammx0 = gmx
            self.betay0 = bay
            self.alfay0 = aly
            self.gammy0 = gmy
            
            # Probe: twiss-functions durch ganze Zelle    
            v_beta=np.array([[bax],[alx],[gmx],[bay],[aly],[gmy]])
            m_cell=self.full_cell.BetaMatrix()
            v_beta_end = m_cell.dot(v_beta)
            print('\nProbe: {Twiss_Ende} == {Zellenmatrix}x{Twiss_Anfang}?')
            print('Anfang: ',v_beta.T)
            print('Ende  : ',v_beta_end.T,'\n')
        
        return self.full_cell,self.betax0,self.betay0
    #-----------------------
    def reverse(self):  ## return a reversed Lattice
        res=Lattice()
        seq=copy(self.seq)
        seq.reverse()
        for ipos in seq:
            elm,s,s=ipos
            res.add_element(elm)
        return res
    #-----------------------
    def append(self,piece):  ## append a Lattice piece
        seq=copy(piece.seq)  
        for ipos in seq:
            elm,s0,s1=ipos
            self.add_element(elm)
    #-----------------------
    def beta_functions(self,steps=10):
        beta_fun=[]
        bx = self.betax0
        ax = self.alfax0
        gx = self.gammx0
        by = self.betay0
        ay = self.alfay0
        gy = self.gammy0
        v_beta=np.array([[bx],[ax],[gx],[by],[ay],[gy]])
        s=0.0
        for ipos in self.seq:
            element,s0,s1 = ipos
            for i_element in element.step_through(steps):
                m_beta = i_element.BetaMatrix()
                v_beta = m_beta.dot(v_beta)
                s += i_element.length
                betax = v_beta[0,0]
                betay = v_beta[3,0]
                viseo = i_element.viseo
                beta_fun.append((s,betax,betay,viseo))
        return beta_fun
        
    def dispersion(self,steps=10,closed=True):
        traj=[]
        v_0=np.array([[0.],[0.],[0.],[0.],[1.]])
        if closed == True:
            m_cell = self.full_cell
            m11=m_cell.matrix[0,0]
            m15=m_cell.matrix[0,4]
            d0 = m15/(1.-m11)     # from H.Wiedemann (6.79) pp.206
            v_0=np.array([[d0],[0.],[0.],[0.],[1.]])
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
###################################################
def test1():
    lattice=SETUP.make_lattice()
    mcell,betax,betay=lattice.cell()
    beta_matrix = mcell.BetaMatrix()
    
    eigen, vectors = LA.eig(beta_matrix)
    print('eigen\n',eigen)
    print('vectors\n',vectors)
    print('Mit numpy.linalg berechneter Eigenwert: \n',eigen[0].real)
    bx=vectors[0][0].real; ax=vectors[1][0].real; gx=vectors[2][0].real
    print('...und sein Eigenvektor dazu: \n {:.6f} {:.6f} {:.6f}'.
          format(bx,ax,gx))
    probe=beta_matrix.dot(vectors[:,0])
    print('Probe: \n {:.6f} {:.6f} {:.6f}'.
          format(probe[0].real,probe[1].real,probe[2].real))
    print(
'''Ergebnis scheint trivial zu sein!
Eigenvvektor ist auf 1 normiert.
Jedes Vielfache eines Eigenvektors ist auch Eigenvektor!
Wo ist der Wert von betax geblieben?
Antwort: Da die Emittanz noch nicht bekannt ist kann man den
Eigenvektor auf 1 normieren und erhält dasselbe wie mit Wille's Formeln.'''
)
    print('==================================================')
    return

def test2():
    lattice=SETUP.make_lattice()
    ## cell boundaries
    mcell,betax,betay=lattice.cell()
    print('BETAx[0] {:.3f} BETAy[0] {:.3f}'.format(betax,betay))
    ## lattice function as f(s)
    beta_fun = lattice.beta_functions(steps=100)   
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
    show()
####################################################
if __name__ == '__main__':
    test1()
    test2()
