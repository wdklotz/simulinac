"""
    Nulstellen durch Intervallschachtelung
    hier: 
        Vergleich Naeherung mit exakter Rechnung fuer stabile Phasen der long. Separatrix
"""
from math import degrees,radians,sin,cos,fabs,pi

def fun(x,phis):
    return (phis+x)*cos(phis)-sin(phis)-sin(x)   # T.Wangler (6.22)
def funp(x,phis):
    return cos(phis)-cos(x)

def intervall(phis):
    a = pi/2.
    b = -a
    # finde Intervall
    while True:
        fa = fun(a,phis)
        fb = fun(b,phis)
        # print((a,b),(fa,fb))
        if (fa>0. and fb<0.) or (fa<0. and fb>0.):
            break
        else:
            b = b-a
    # finde Nulstelle
    condition = True
    while condition:
        c= (a+b)/2.
        fc = fun(c,phis)
        if (fa>0. and fc>0.) or (fa<0. and fc<0.):
            a = c
            fa = fun(a,phis)
        if (fb>0. and fc>0.) or (fb<0. and fc<0.):
            b = c
            fb = fun(b,phis)
        # print((a,b),(fa,fb))
        if fabs(fa-fb) < 1.e-7:
            condition = False
    return (a+b)/2.
def newton(phis):
    condition = True
    phin = phis-1.e-3
    while condition:
        fn     = fun(phin,phis)
        phinp1 = phin - fn/funp(phin,phis)
        # print((phin,phinp1,fun(phinp1,phis)))
        if fabs(fun(phinp1,phis))<1.e-7:
            condition = False
        else:
            phin = phinp1
    return phinp1
    
# Ergebnis
phis = radians(-30.)
print('Fuer synchrone Phase {} [deg]'.format(degrees(phis)))
print('---------Intervall----------')
phi2 = intervall(phis)
print('intervall:  phi2 {} psis {} [deg]'.format(degrees(phi2),degrees(fabs(phi2)+fabs(phis))))
print('---------Newton-------------')
phi2 = newton(phis)
print('newton:     phi2 {} psis {} [deg]'.format(degrees(phi2),degrees(fabs(phi2)+fabs(phis))))
print('T.Wrangler: phi2 {} psis {} [deg] (pp.178)'.format(degrees(2.*phis),degrees(3.*fabs(phis))))

