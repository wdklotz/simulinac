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
from elements import D,QFth,QDth,QF,QD
from lattice import Lattice
from numpy import linalg as LA
from matplotlib.pyplot import plot,show,legend,figure,subplot,axis
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
    show(block=False)
def make_thin (kf1,kf2,ld,anz=1,verbose=False):
    kd1 = kf1       # k1 for qf1 & qd1
    ld1 = lf1 = 0.4 # len qf1 & qd1
    ff1 = kf1*lf1   # focal length qf1
    fd1 = kd1*ld1   # focal length qd1

    kf2 = kf2
    kd2 = kf2       # same as above for qf2 & qd2
    lf2 = ld2 = lf1
    ff2 = kf2*lf2
    fd2 = kd2*ld2

    slices = 1      # present the quad by 1 thin-quads   (bad!)
    slices = 3      # present the quad by 3 thin-quads   (better!)
    slices = 6      # present the quad by 6 thin-quads   (near perfect!)
    
    DL  = D(length=ld,label='L')
    QF1 = QFth(k0=kf1,length=0.5*lf1/slices,label='QF1')
    QF2 = QFth(k0=kf2,length=0.5*lf2/slices,label='QF2')
    QD1 = QDth(k0=kd1,length=0.5*ld1/slices,label='QD1')
    QD2 = QDth(k0=kd2,length=0.5*ld2/slices,label='QD2')

    cell = Lattice()
    for i in range(slices): cell.add_element(QD1)
    for i in range(slices): cell.add_element(QF1)
    cell.add_element(DL)
    for i in range(slices): cell.add_element(QD2)
    for i in range(slices): cell.add_element(QF2)
    for i in range(slices): cell.add_element(QF2)
    for i in range(slices): cell.add_element(QD2)
    cell.add_element(DL)
    for i in range(slices): cell.add_element(QF1)
    for i in range(slices): cell.add_element(QD1)

    lat = Lattice()
    for i in range(anz):
        lat.concat(cell)
    mcell,betax,betay=lat.cell()

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
        lat.concat(cell)
    mcell,betax,betay=lat.cell()
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
    print('---------------------------------Test0---')
    found=search()
    print('found minimal with: (kf, kd, L)= ',found)
def test1(kf,kd,ld):
    print('---------------------------------Test1---')
    print('using kf,kd,ld',kf,kd,ld)
    cell,dummy,dummy = make_thin(kf,kd,ld)
    mcell,betax,betay=cell.cell()
    beta_matrix = mcell.beta_matrix()
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
    print('---------------------------------Test2---')
    print('using kf,kd,ld',kf,kd,ld)
    cell,dummy,dummy = make_thick(kf,kd,ld)
    mcell,betax,betay=cell.cell()
    beta_matrix = mcell.beta_matrix()
    eigen, vectors = LA.eig(beta_matrix)
    print('eigen\n',eigen)
    print('vectors\n',vectors)
    print('Mit numpy.linalg berechneter Eigenwert: \n',eigen[0].real)
    bx=vectors[0][0].real; ax=vectors[1][0].real; gx=vectors[2][0].real
    print('...und sein Eigenvektor dazu: \n {:.6f} {:.6f} {:.6f}'.format(bx,ax,gx))
    probe=beta_matrix.dot(vectors[:,0])
    print('Probe: \n {:.6f} {:.6f} {:.6f}'.format(probe[0].real,probe[1].real,probe[2].real))
    print("""
    ================================
    I don't understand this result!!
    ================================
    """)
    cell.symplecticity()
def test3(kf,kd,ld):
    print('---------------------------------Test3---')
    print('using kf,kd,ld',kf,kd,ld)
    anz = 3
    # thin
    cell,dummy,dummy    = make_thin(kf,kd,ld,anz=anz)
    mcell,betax,betay   = cell.cell()
    beta_matrix         = mcell.beta_matrix()
    beta_fun_thin       = cell.twiss_functions(steps=100)
    # thick
    cell,dummy,dummy    = make_thick(kf,kd,ld,anz=anz)
    mcell,betax,betay   = cell.cell()
    beta_matrix         = mcell.beta_matrix()
    beta_fun_thick      = cell.twiss_functions(steps=100)
    display((beta_fun_thin,beta_fun_thick))
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    # test0()
    # test1(5.,5.,1.)
    # test2(5.,5.,1.)
    # test3(5.,5.,1.)     # gesund!
    test3(4.,4.,1.2)    # gesünder!
    # test3(6.,6.,1.8)      # extrem!
