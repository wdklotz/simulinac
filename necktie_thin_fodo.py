#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# -*- coding: utf-8 -*-
"""
Copyright 2016 Wolf-Dieter Klotz <wdklotz@gmail.com>
"""
from math import sqrt, sin, cos, sinh, cosh
import matplotlib.pyplot as plt

def necktie(u,v):
    """
        function checks if (u,v) lies in stable necktie for asymmetric FODO in thin lens approximation. (see Wiedemann pp. 193 Vol 1)
        u = L/f1
        v = L/f2
        returns True if stable
        note: if u>0 then v<0 and vice versa for a thin lens FODO
    """
    wx=u+v-u*v
    wy=-u-v-u*v # change signs of u,v for 2nd tranverse plane
    if (0<wx and wx<1) and (0<wy and wy<1):
        return True
    else:
        return False

def betagamma(tkin):    # beta*gamma for proton as function of kin. energy
    tkin_zu_e0 = tkin/938.  # tkin in [MeV]
    res = 2.*tkin_zu_e0 + tkin_zu_e0**2
    res = sqrt(res)
    return res

def kq(grad,tkin):   # quad k [1/m^2] for proton with beta*gamma momentum and grad[T/m] gradient
    res = 0.31952*grad/betagamma(tkin)
    return res

##-------------------------------------Test0--
def test0():
    uwerte=[0.+i*0.01 for i in range(int(1./0.01))]
    print(uwerte)
    vwerte=[-i for i in uwerte]

    stabile_werte=[]
    for u in uwerte:
        for v in vwerte:
            if necktie(u,v): stabile_werte.append((u,v))

    stabile_uwerte=[abs(t[0]) for t in stabile_werte]
    stabile_vwerte=[abs(t[1]) for t in stabile_werte]
    plt.scatter(stabile_uwerte,stabile_vwerte)
    plt.show(block=True)

##-------------------------------------Test1--
def test1():
    """ quad gradient scan """
    ld=1.
    lq=0.1
    L=ld+lq
    tkin = 1.
    grad1 = +20.
#     grad2 = -20.
    ssize = 0.04
    bg = betagamma(tkin)

    grad1_werte=[0.1+i*ssize for i in range(int(grad1/ssize))]
    grad2_werte=[-i for i in grad1_werte]
    stabile_werte=[]
    for g1 in grad1_werte:
        u=L*kq(g1,tkin)*lq
        for g2 in grad2_werte:
            v=L*kq(g2,tkin)*lq
            if necktie(u,v):
                stabile_werte.append((g1,g2))

    stabile_g1_werte = [abs(t[0]) for t in stabile_werte]
    stabile_g2_werte = [abs(t[1]) for t in stabile_werte]
    plt.scatter(stabile_g1_werte,stabile_g2_werte)
    plt.title("axes are dB/dx[T/m], Tk={:4.4}[MeV], L={:4.4}[m]".format(tkin,L))
    plt.show(block=True)


if __name__ == '__main__':
     test0()
#test1()