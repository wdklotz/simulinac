#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import dictprnt,objprnt
from elements import D,QFth,QDth,QF,QD
from lattice import Lattice
from numpy import linalg as LA
from pylab import plot,show,legend,figure,subplot,axis
from math import fabs

def display(functions):
    bthin = functions[0]
    bthick = functions[1]
    s  = [x[0] for x in bthin]    # s
    xs = [x[1] for x in bthin]    # betax
    ys = [x[2] for x in bthin]    # betay
    viseo = [x[3] for x in bthin]
    zero  = [0. for x in bthin]  # zero line
    plot(s,xs,label='bx/thin')
    plot(s,ys,label='by/thin')
    
    s  = [x[0] for x in bthick]    # s
    xs = [x[1] for x in bthick]    # betax
    ys = [x[2] for x in bthick]    # betay
    viseo = [x[3] for x in bthick]
    zero  = [0. for x in bthick]  # zero line
    plot(s,xs,label='bx')
    plot(s,ys,label='by')

    vscale=axis()[3]*0.1
    viseo = [x*vscale for x in viseo]
    plot(s,viseo,label='element',color='black')
    plot(s,zero,color='black')
    legend(loc='lower right',fontsize='x-small')
    show(block=True)
def make_thin (kf1,kf2,ld,anz=1,verbose=False):
    kf1 = kf1
    kd1 = kf1
    lf1 = 0.4
    ld1 = lf1
    ff1 = kf1*lf1
    fd1 = kd1*ld1

    kf2 = kf2
    kd2 = kf2
    lf2 = lf1
    ld2 = lf2
    ff2 = kf2*lf2
    fd2 = kd2*ld2
            
    ld = ld
    
    DL  = D(length=ld,label='L')
    QF1 = QFth(k0=kf1,length=0.5*lf1,label='QF1')
    QF2 = QFth(k0=kf2,length=0.5*lf2,label='QF2')
    QD1 = QDth(k0=kd1,length=0.5*ld1,label='QD1')
    QD2 = QDth(k0=kd2,length=0.5*ld2,label='QD2')
    
    cell = Lattice()
    cell.add_element(QD1)
    cell.add_element(QF1)
    cell.add_element(DL)
    cell.add_element(QD2)
    cell.add_element(QF2)
    cell.add_element(QF2)
    cell.add_element(QD2)
    cell.add_element(DL)
    cell.add_element(QF1)
    cell.add_element(QD1)
    lat = Lattice()
    for i in range(anz):
        lat.append(cell)
    mcell,betax,betay=lat.cell(verbose=verbose)
    if verbose:
        # {:.3f}
        print('L= {:.3f}'.format(ld),end=' ')
        print('triplet 1: kf= {:.3f}, kd={:.3f}'.format(kf1,kd1),end=' | ')
        print('triplet 2: kf= {:.3f}, kd={:.3f}'.format(kf2,kd2),end=' | ')
        print('betax= {:.3f}, betay= {:.3f}'.format(betax,betay))
    return lat,betax,betay
def make_thick(kf1,kf2,ld,anz=1,verbose=False):
    kf1 = kf1
    kd1 = kf1
    lf1 = 0.4
    ld1 = lf1
    ff1 = kf1*lf1
    fd1 = kd1*ld1

    kf2 = kf2
    kd2 = kf2
    lf2 = lf1
    ld2 = lf2
    ff2 = kf2*lf2
    fd2 = kd2*ld2
            
    ld = ld
    
    DL  = D(length=ld,label='L')
    QF1 = QF(k0=kf1,length=0.5*lf1,label='QF1')
    QF2 = QF(k0=kf2,length=0.5*lf2,label='QF2')
    QD1 = QD(k0=kd1,length=0.5*ld1,label='QD1')
    QD2 = QD(k0=kd2,length=0.5*ld2,label='QD2')
    
    cell = Lattice()
    cell.add_element(QD1)
    cell.add_element(QF1)
    cell.add_element(DL)
    cell.add_element(QD2)
    cell.add_element(QF2)
    cell.add_element(QF2)
    cell.add_element(QD2)
    cell.add_element(DL)
    cell.add_element(QF1)
    cell.add_element(QD1)
    lat = Lattice()
    for i in range(anz):
        lat.append(cell)
    mcell,betax,betay=lat.cell(verbose=verbose)
    if verbose:
        # {:.3f}
        print('L= {:.3f}'.format(ld),end=' ')
        print('triplet 1: kf= {:.3f}, kd={:.3f}'.format(kf1,kd1),end=' | ')
        print('triplet 2: kf= {:.3f}, kd={:.3f}'.format(kf2,kd2),end=' | ')
        print('betax= {:.3f}, betay= {:.3f}'.format(betax,betay))
    return lat,betax,betay
def search():
    xmin = ymin = 1.e25
    xmax = ymax = 0.
    crit = xmin
    for kf in [4.0+n*0.5 for n in range(10)]:
        for kd in [4.0+n*0.5 for n in range(10)]:
            for ld in [1.0+n*0.2 for n in range(10)]:
                try:
                    m,bx,by = make_thin(kf,kd,ld)
                except RuntimeError:
                    continue
                else:
                    # x=bx*by
                    if bx < xmin:
                        xmin=bx
                    if by < ymin:
                        ymin = by
                    if bx > xmax:
                        xmax = bx
                    if by > ymax:
                        ymax = by
                    x=fabs(xmax-xmin)+fabs(ymax-ymin)
                    if x < crit and x != 0.:
                        print(x)
                        crit = x
                        found=(kf,kd,ld) 
    return found 
def test0():
    found=search()
    print('found minimal with: (kf, kd, L)= ',found)
def test1(kf,kd,ld):
    print('test1: using kf,kd,ld',kf,kd,ld)
    cell,dummy,dummy = make_thin(kf,kd,ld)
    mcell,betax,betay=cell.cell(verbose=True)
    beta_matrix = mcell.BetaMatrix()    
    eigen, vectors = LA.eig(beta_matrix)
    print('eigen\n',eigen)
    print('vectors\n',vectors)
    print('Mit numpy.linalg berechneter Eigenwert: \n',eigen[0].real)
    bx=vectors[0][0].real; ax=vectors[1][0].real; gx=vectors[2][0].real
    print('...und sein Eigenvektor dazu: \n {:.6f} {:.6f} {:.6f}'.format(bx,ax,gx))
    probe=beta_matrix.dot(vectors[:,0])
    print('Probe: \n {:.6f} {:.6f} {:.6f}'.format(probe[0].real,probe[1].real,probe[2].real))
    cell.symplecticity()
def test2(kf,kd,ld):
    print('test2: using kf,kd,ld',kf,kd,ld)
    cell,dummy,dummy = make_thick(kf,kd,ld)
    mcell,betax,betay=cell.cell(verbose=True)
    beta_matrix = mcell.BetaMatrix()    
    eigen, vectors = LA.eig(beta_matrix)
    print('eigen\n',eigen)
    print('vectors\n',vectors)
    print('Mit numpy.linalg berechneter Eigenwert: \n',eigen[0].real)
    bx=vectors[0][0].real; ax=vectors[1][0].real; gx=vectors[2][0].real
    print('...und sein Eigenvektor dazu: \n {:.6f} {:.6f} {:.6f}'.format(bx,ax,gx))
    probe=beta_matrix.dot(vectors[:,0])
    print('Probe: \n {:.6f} {:.6f} {:.6f}'.format(probe[0].real,probe[1].real,probe[2].real))
    cell.symplecticity()
def test3(kf,kd,ld):
    print('test3: using kf,kd,ld',kf,kd,ld)
    anz = 3
    # thin 
    cell,dummy,dummy = make_thin(kf,kd,ld,anz=anz)
    mcell,betax,betay=cell.cell(verbose=True)
    beta_matrix = mcell.BetaMatrix()    
    beta_fun_thin,cl,sl = cell.functions(steps=100)   
    # thick
    cell,dummy,dummy = make_thick(kf,kd,ld,anz=anz)
    mcell,betax,betay=cell.cell(verbose=True)
    beta_matrix = mcell.BetaMatrix()    
    beta_fun_thick,cl,sl = cell.functions(steps=100)   
    display((beta_fun_thin,beta_fun_thick))
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    # test0()
    # test1(5.,5.,1.)
    # test2(5.,5.,1.)
    # test3(5.,5.,1.)     # gesund!
    test3(4.,4.,1.2)     # gesünder!
    # test3(7.,7.,2.1)  # extrem!
