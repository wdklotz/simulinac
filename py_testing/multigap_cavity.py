import matplotlib.pyplot as PLT
import math as M
import numpy as NP
from scipy.stats import norm
import scipy.constants as C

pi      = M.pi
twopi   = pi+pi
# threepi = twopi+pi
# fourpi  = threepi+pi
 
# cl    = C.c
# freq  = 816e6
# lamb  = cl/freq
# omega = twopi*freq
# m0c2  = C.value('proton mass energy equivalent in MeV')
# tk    = 50   # kin. energy
# gamma =  1.+tk/m0c2
# beta  = M.sqrt(1. - 1./gamma**2)
# blamb = beta*lamb

# print(" \
#     c-light [m/s]....... {0:G} \n \
#     freq [Hz]..........  {1:G} \n \
#     omega [Hz].........  {3:G} \n \
#     rest mass [MeV]..... {4:G} \n \
#     kin energy [Mev].... {5:G} \n \
#     gamma []............ {6:G} \n \
#     beta []............. {7:G} \n \
#     lambda [m].......... {2:G} \n \
#     beta-lambda [m]..... {8:G} \
# ".format(cl,freq,lamb,omega,m0c2,tk,gamma,beta,blamb))
# print(omega*lamb/cl/twopi)

def generate_space_part(delta=0.0):
    gn,gsig     = (5,  0.125)   # nb gaps, gap sigma
    glen        = 1.0           # interval length
    gstart      = -glen/2.
    # concatenate several gaps
    y = NP.array([])
    x = NP.array([])
    xstart = gstart
    for n in range(gn):
        xstop  = xstart + glen*(1.+delta/100.)
        xrange = NP.linspace(xstart, xstop, 50)
        xpos   = (xstart+glen/2.)
        yy     = M.sqrt(twopi)*gsig * norm.pdf(xrange,xpos,gsig)
        x      = NP.concatenate((x,xrange))
        y      = NP.concatenate((y,yy))        # space part
        xstart = xstop
    return x,y

phi = lambda u: twopi*u+p0

fig = PLT.figure(figsize=(14.,6.))      # prep plotting

phisoll = [-45.0,-90.0,-0.0]
x,y = generate_space_part(delta=0.0)    # gap positions
splt = PLT.subplot(211)
PLT.title(r'5 gaps with spacing = $\beta\lambda$')
PLT.plot(x,y,':')
for p in phisoll:
    p0 = M.radians(p)
    Eup = NP.array([])
    for i,u in enumerate(x):
        v = M.cos(phi(u))*y[i]
        Eup = NP.append(Eup,[v])
    fmt = 'ok' if p == phisoll[0] else '-'
    PLT.plot(x,Eup,fmt,label=f'{p:+3.1f} deg',markersize=3)
PLT.legend()

x,y = generate_space_part(delta=+10.0)    # gap positions +10%
splt = PLT.subplot(212)
PLT.title(r'5 gaps with spacing = $\beta\lambda$ + 10%')
PLT.plot(x,y,':')
for p in phisoll:
    p0 = M.radians(p)
    Eup = NP.array([])
    for i,u in enumerate(x):
        v = M.cos(phi(u))*y[i]
        Eup = NP.append(Eup,[v])
    fmt = 'ok' if p == phisoll[0] else '-'
    PLT.plot(x,Eup,fmt,label=f'{p:+3.1f} deg',markersize=3)
PLT.legend()
splt.set_xlabel(r'u=s/$\beta\lambda$')

PLT.show()