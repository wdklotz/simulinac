#!/Users/klotz/pyzo2015a/python
# -*- coding: utf-8 -*-
from setup import Phys,Proton
from pylab import plot,show,legend
from math import cos,pi,sqrt
'''try Ruth's symplectic integrator'''
def display(werte):
    i = [x[0] for x in werte]
    p = [x[1] for x in werte]
    q = [x[2] for x in werte]
    t = [x[3] for x in werte]
    plot(t,p,label='p')
    plot(t,q,label='q')
    # plot(q,p,label='w(phi)')
    legend(loc='lower right',fontsize='x-small')
    show()
def step1(p0,q0,s0):
    p1 = p0
    q1 = q0 + 0.5 * h * p1
    s1 = s0 + h
    return (p1,q1,s1)
def step2(p1,q1,s1):
    t = s1 
    p = p1 + t * B *(cos(q1)+cos(qs))
    q = q1 + 0.5 * h * p
    return (p,q,t)
#-----------*-----------*-----------*-----------*-----------*-----------*-----------*
rad=Phys['radians']
deg=Phys['degrees']
Ws=50.            # soll energie
prot=Proton(Ws)
bs=prot.beta
gs=prot.gamma
Ts=prot.TrTf()
m0=prot.e0        # ruhmasse
E0=Phys['spalt_spannung']/Phys['spalt_laenge']       # feldstärke
# Ez=Phys['spalt_spannung']                          # feldstärke
lamb=Phys['wellenlänge']        # wellenlänge
A=2.*pi/(bs*gs)**3/lamb   # faktor A
B=E0*Ts/m0   # faktor B
wA=sqrt(A)      # wurzel A
phis=-45.   # soll phase
phi=phis-5.
print('phis\t{:.4g}'.format(phis))
print('phi\t{:.4g}'.format(phi))
phi =rad*phi  # akt. phase
phis=rad*phis
w=-0.0/(m0/Ws)         # akt. delta-w [%]/mc2 [%]

p0 = wA*w    # p kanonisch
q0 = phi  # q kanonish
qs = phis
s0 = 0.   # s==t als unabhängige variable
h  = 0.02 # step size of integrator

print('E0\t{:.4g}'.format(E0))
print('A\t{:.4g}'.format(A))
print('B\t{:.4g}'.format(B))
print('wA\t{:.4g}'.format(wA))
print('B/wA\t{:.4g}'.format(B/wA))
print('p0\t{:.4g}'.format(p0))
print('q0\t{:.4g}'.format(q0))
print('qs\t{:.4g}'.format(qs))
print('h\t{:.4g}'.format(h))

print('-----------*-------*-----------')
print('unverständliches Ergenis hier!!')
print('-----------*-------*-----------')

werte=[]
for i in range(10000):
    (p1,q1,s1) = step1(p0,q0,s0)
    (p ,q ,s ) = step2(p1,q1,s1)
    werte.append([i,p/wA,(q-qs)*deg,s])
    p0 = p
    q0 = q
    s0 = s
    # print(i,p,q,s)
display(werte)

    
    
    
