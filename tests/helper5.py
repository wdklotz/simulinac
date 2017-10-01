import sys
sys.path.insert(0,'..')

from math import sqrt
from setutil import tblprnt

c     = 2.997e8  #m/s
m0    = 938.26   #MeV
freq  = 815.0e6  #Hz
lamb  = c/freq   #m
z     = 0.044    #m

def delta_phi(T,z):
    gamma = (T+m0)/m0
    beta  = sqrt(1.-1./gamma**2)
    phi   = 360/beta/lamb*z
    return phi

header = ['T[Mev]','z[m]','delta phi[deg]','delta phi[deg/cm]']
tablrows = []
for t in [10.,50.,100.,150.,200.,250.,1.e3]:
    dphi = delta_phi(t,z)
    tablrow = ['{:8.4g}'.format(t),
                '{:8.4f}'.format(z),
                '{:8.4f}'.format(dphi),
                '{:8.4f}'.format(dphi/(z*1.e2))]
    tablrows.append(tablrow)
print('Phase advance as function of kinetic energy T and distance z')
print('c={:8.4g}[m/s], m0={:8.4g}[MeV/c2], fr={:8.4g}[MHz]'.format(c,m0,freq*1.e-6))
print(tblprnt(header,tablrows))