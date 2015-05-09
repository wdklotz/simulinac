#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import Phys,Proton
from pylab import plot,show,legend,figure,subplot,axis
from math import pi,cos

def display(functions):
    for func in functions:
        offp=  [x[0] for x in func]
        offw=  [x[1] for x in func]
        plot(offp,offw,label=r'$\Delta$w($\Delta$phi)')
    # legend(loc='lower right',fontsize='x-small')
    show()
def phase(factor,ws,w):
    return  (w-ws)*factor
def energy(factor,ps,p):
    return factor*(cos(ps)-cos(p))
def loop(phis0,ws0,phi0,w0):
    ps=phis0
    p=phi0
    ws=ws0
    w=w0
    
    deg = Phys['degrees']
    pr=Proton(ws); mc2=pr.e0; g=pr.gamma; b=pr.beta; phase_factor = -2*pi/(mc2*g*g*g*b*b)
    T=pr.TrTf(); energy_factor=Phys['spalt_spannung']*T

    func=[]
    for i in range(70):
        # print('{}   (ps,ws)= ({:5.3f},{:5.3f})  (p,w)= ({:5.3f},{:5.3f})'.format(i,ps*deg,ws,p*deg,w))
        w = w + energy(energy_factor,ps,p)
        p = p + phase(phase_factor,ws,w)
        fun=((ps-p)*deg,w-ws)
        func.append(fun)
    return func
def test0():
    rad = Phys['radians']
    phis0=-30.
    ws0=50.
    w0 =ws0
    functions=[]
    for i in [phis0+n*4.775065 for n in range(20)]:
        phi0 = i*rad
        func= loop(phis0*rad,ws0,phi0,w0)
        functions.append(func)
    display(functions)
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    test0()
        
        