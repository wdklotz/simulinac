# -*- coding: utf-8 -*-
from setup import fodo
from elements import I,D,QF,QD,SD,WD
from lattice import Lattice

from pylab import plot, show, legend

def plotter(beta_fun):
    s  = [x[0] for x in beta_fun]
    xs = [x[1] for x in beta_fun]
    ys = [x[2] for x in beta_fun]
    vs = [x[3]-1. for x in beta_fun]
    zero=[-1. for x in beta_fun]
##    for i in range(len(s)):
##        print('s, betax(s) betay(s)',s[i],xs[i],ys[i])
    plot(s,xs,label='betax')
    plot(s,ys,label='betay')
    plot(s,vs,label='element',color='black')
    plot(s,zero,color='black')
    legend(loc='upper left')
    show()

def test1():
    print("cloned from K.Wille's Beispiel auf pp. 112-113 (wdk)")
    kqf=  fodo()['k_quad_f']
    lqf=  fodo()['length_quad_f']
    kqd=  fodo()['k_quad_d']
    lqd=  fodo()['length_quad_d']
    rhob= fodo()['beding_radius']
    lb=   fodo()['dipole_length']
    ld=   fodo()['drift_length']
    ## elements
    mqf=QF(kqf,lqf,'QF')
    mqd=QD(kqd,lqd,'QD')
    mb=SD(rhob,lb,'B')
    mw=WD(lb,rhob)
    md=D(ld)    
    ## lattice
    lattice=Lattice()
    lattice.add_element(mqf)
    lattice.add_element(md)
    lattice.add_element(mw)
    lattice.add_element(mb)
    lattice.add_element(mw)
    lattice.add_element(md)
    lattice.add_element(mqd)
    lattice.add_element(md)
    lattice.add_element(mw)
    lattice.add_element(mb)
    lattice.add_element(mw)
    lattice.add_element(md)
    lattice.add_element(mqf)
    lattice.out()
    ## cell boundaries
    mcell,betax,betay=lattice.cell()
    print('BETAx[0] {:.3f} BETAy[0] {:.3f}'.format(betax,betay))
    ## lattice function as f(s)
    beta_fun = lattice.beta_functions(120)    
    ## plots
    plotter(beta_fun)
####################################################
if __name__ == '__main__':
    test1()
    
