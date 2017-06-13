#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# -*- coding: utf-8 -*-
"""
Copyright 2016 Wolf-Dieter Klotz <wdklotz@gmail.com>
"""
import sys
from math import sqrt, sin, cos, sinh, cosh
from matplotlib.pyplot import scatter,show,legend,figure,subplot,axis

def RM44_thin(ff,fd,L):
    f1 = ff
    f2 = -fd
    a = 1./f1+1./f2-L/(f1*f2)
    fstar = 1./a
    R=[[0 for i in range(4)] for j in range(4)]
    R[0][0] = 1.-2.*L/fstar
    R[0][1] = 2.*L*(1.-L/f2)
    R[1][0] = -2./fstar*(1.-L/f1)
    R[1][1] = 1.-2.*L/fstar
    R[2][2] = 1.-2.*L/fstar
    R[2][3] = 2.*L*(1.-L/f1)
    R[3][2] = -2./fstar*(1.-L/f2)
    R[3][3] = 1.-2.*L/fstar
    Trxx = R[0][0]+R[1][1]
    Tryy = R[2][2]+R[3][3]
    return (R,Trxx,Tryy)

def RM44(kd,kf,ld,lq):
    """
    kd: D-quad k value in [1/mˆ2]
    kf: F-quad k value in [1/mˆ2]
    ld: drift lentbth between quads in [m]
    lq: quad length in [m]

    returns a tuple (R,Trxx,Tryy)
        R: completw lin matrix of a FODO cell
        Trxx: trace of Rxx submatrix
        Tryy: trace of Ryy submatrix
    """
    global params
    nd = 4
    R =[[0 for i in range(nd)] for j in range(nd)]
    QD=[[0 for i in range(nd)] for j in range(nd)]
    QF=[[0 for i in range(nd)] for j in range(nd)]
    DR=[[0 for i in range(nd)] for j in range(nd)]

    klq  = sqrt(kd)*lq         # [dimensionless]
    QD[0][0] = cosh(klq)
    QD[0][1] = sinh(klq)/klq
    QD[1][0] = klq*sinh(klq)
    QD[1][1] = cosh(klq)
    QD[2][2] = cos(klq)
    QD[2][3] = sin(klq)/klq
    QD[3][2] = -klq*sin(klq)
    QD[3][3] = cos(klq)

#     out('QD',(QD,0.,0.))

    klq  = sqrt(kf)*lq         # [dimensionless]
    QF[0][0] = cos(klq)
    QF[0][1] = sin(klq)/klq
    QF[1][0] = -klq*sin(klq)
    QF[1][1] = cos(klq)
    QF[2][2] = cosh(klq)
    QF[2][3] = sinh(klq)/klq
    QF[3][2] = klq*sinh(klq)
    QF[3][3] = cosh(klq)

#     out('QF',(QF,0.,0.))

    DR[0][0] = 1.
    DR[0][1] = ld
    DR[1][1] = 1.
    DR[2][2] = 1.
    DR[2][3] = ld
    DR[3][3] = 1.

#     out('DR',(DR,0.,0.))

    latt = [DR,QD,DR,QF]    # lattice
    R = latt[0]
    for i in range(1,len(latt)):
        R = mult(latt[i],R)

    Trxx = R[0][0]+R[1][1]
    Tryy = R[2][2]+R[3][3]
    return (R,Trxx,Tryy)

def mult(A,B):          # multplies marices A*B
    ra = len(A[0])
    ca = len(A[1])
    rb = len(B[1])
    cb = len(B[0])
    R = [[0 for i in range(ra)] for j in range(cb)]
    for i in range(ra):
     for j in range(cb):
        for k in range(ca):
            R[i][j] += A[i][k]*B[k][j]
    return R

def out(what,matrix):       # print a matrix with title
    R = matrix[0]
    Trxx = matrix[1]
    Tryy = matrix[2]
    print('{:s}   cos(mux)={:8.4f} cos(muy)={:8.4}'.format(what,abs(Trxx)/2.,abs(Tryy)/2.))
    for i in range(len(R[0])):
        for j in range(len(R[1])):
            if j < len(R[1])-1:
                print('{:8.4f}  '.format(R[i][j]),end='')
            else:
                print('{:8.4f}  '.format(R[i][j]))

def betagamma(tkin):    # beta*gamma for proton as function of kin. energy
    tdive0 = tkin/938.  # tkin in [MeV]
    res = 2.*tdive0 + tdive0**2
    res = sqrt(res)
    return res

def kq(g,tkin):   # quad k [1/m^2] for proton with beta*gamma momentum and g [T/m] gradient
    res = 0.31952*g/betagamma(tkin)
    return res

params = {}
stabile_werte=[]        # list stable solutions

ld = 1.     # l drift [m]
lq = 0.1    # l quad [m]
tkin = 5.  # T kin [MeV]
L = ld+lq
gdmx = 20.    # B gradient [T/m]
gfmx = 20.    # B gradient [T/m]
bgamma = betagamma(tkin)

stepsize = 0.05
print(range(int(gdmx/stepsize)))
for gd in [0.01+i*stepsize for i in range(int(gdmx/stepsize))]:     #vary D=quad gradient
    for gf in [0.01+j*stepsize for j in range(int(gfmx/stepsize))]: #vary F-quad gradient
        kd = kq(gd,tkin)    # [1./m^2]
        kf = kq(gf,tkin)    # [1./m^2]
        fd = 1./(kd*lq)     # thin lens foc length
        ff = 1./(kf*lq)     # thin lens foc length
        u = L/fd
        v = L/ff

#         matrix = RM44(kd,kf,ld,lq)
        matrix = RM44_thin(ff,fd,L)
        Trxx = matrix[1]
        Tryy = matrix[2]
        cosmux = abs(Trxx)/2.
        cosmuy = abs(Tryy)/2.
        if cosmux > 1. and cosmuy > 1.: continue       # discard unstable

        stabile_werte.append((gd,kd,cosmux,gf,kf,cosmuy,u,v))

#         print('==================================')
#         params['drift btw. quads[m]'] = ld
#         params['quad len[m]'] = lq
#         params['thin lens f[m]'] = fx
#         params['kinetic T[Mev]'] = tkin
#         params['QD gradient g[T/m]'] = gd
#         params['QF gradient g[T/m]'] = gf
#         params['QD strenth k[1/m**2]'] = kd
#         params['QF strenth k[1/m**2]'] = kf
#         params['beta*gamma'] = bgamma
#         for i in sorted(params.items()):
#             print ('|{:s}= {:4.4}| '.format(i[0],i[1]),end='')
#         print('')

#         out('thick',matrix)
#         print('f-foc={:4.4} f-defoc={:4.4}'.format(ff,fd))
#         out('thin',matrix)
# print('stabile Werte')
# for i in stabile_werte:
#     print("gd {:4.4}[T/m] kd {:4.4}[1/mˆ2] cos(mux) {:4.4} | gf {:4.4}[T/m] kf {:4.4}[1/mˆ2] cos(muy) {:4.4}".format(i[0],i[1],i[2],i[3],i[4],i[5]))

def plot_stabile_werte(werte):
    ax_x = [x[6] for x in werte]
    ax_y = [x[7] for x in werte]
    scatter(ax_x,ax_y)
    show(block=True)

plot_stabile_werte(stabile_werte)

