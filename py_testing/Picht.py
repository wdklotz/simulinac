import numpy as NP
import math as MATH

def M(arg):
    a,b,c,d = arg
    return NP.array([[a,b],[c,d]])

def Minverse(arg):
    a,b,c,d = arg
    m = NP.array([[d,-b],[-c,a]])
    det = a*d-b*c
    return m/det

def Mtranspose(arg):
    a,b,c,d = arg
    m = NP.array([[a,c],[b,d]])
    return m
    
arg = [1,0,5,1]
m1 = M(arg)
m2 = Minverse(arg)
m3 = Mtranspose(arg)
print('M\n',m1)
print('Minv\n',m2)
print('Mtrans\n',m3)
print('M*Minv\n',NP.dot(m1,m2))
print('M*Mtrans\n',NP.dot(m1,m3))

def Picht(gamma,inv=False):
    sgn = 1. if inv == False else -1.
    g2m1 = gamma**2-1.
    c=sgn*0.5*gamma/g2m1
    return NP.array([[1,0],[c,1]])*MATH.pow(g2m1,sgn*0.25)
    
T=50.09747033226243
gamma = 1.+T/938.272
pt = Picht(gamma)
pti = Picht(gamma,True)
print('P\n',pt)
print('Pinv\n',pti)
print('P*Pinv\n',NP.dot(pt,pti))

xv=[0.01,0.001]
Xv=NP.dot(pt,xv)
xvr=NP.dot(pti,Xv)
print('x\n',xv)
print('P*x\n',Xv)
print('Pinv*(P*x)\n',xvr)

