#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import Phys,Proton
from pylab import plot,show,legend,figure,subplot,axis
from math import pi,cos

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

    functions=[]
    for i in range(70):
        # print('{}   (ps,ws)= ({:5.3f},{:5.3f})  (p,w)= ({:5.3f},{:5.3f})'.format(i,ps*deg,ws,p*deg,w))
        w = w + energy(energy_factor,ps,p)
        p = p + phase(phase_factor,ws,w)
        functions.append([i,(ps-p)*deg,w-ws])
    bin=   [x[0] for x in functions]
    offp=  [x[1] for x in functions]
    offw=  [x[2] for x in functions]
    plot(offw,offp,label=r'$\Delta$w($\Delta$phi)')
    legend(loc='lower right',fontsize='x-small')
    show()
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
if __name__ == '__main__':
    phis0=Phys['radians']*-90.
    phi0 =Phys['radians']*-95.
    ws0=50.
    w0 =ws0
    loop(phis0,ws0,phi0,w0)
        
        
        