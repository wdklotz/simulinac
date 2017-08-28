from math import sin,cos,radians,pi
from matplotlib.pyplot import plot,show,legend,figure,subplot,axis

"""
Iterate a sin or cos from difference equations
"""
# Folgerung:
#     wenn a12 gleich 0 keine Oszilation
#     wenn a12 ungleich 0 Oszilation; a12 bestimmt anscheinend oder beinflusst die Frequenz
#     a21 muss negativ sein damit Oszilation
def func(yn,h):
    # a11 = a22 = 1.      # fuer reinen sin od. cos
    # a12 = h             # fuer reinen sin od. cos
    l = 0.1               # fuer sumulation T3D
    a11 = a22 = 1.-l*h    # fuer sumulation T3D
    a12 = l*(2.-h*l)      # fuer sumulation T3D !! wenn a12=0. keine oszilation !!
    a21 = -h
    x   = yn[2]
    y   = a11*yn[0] + a12*yn[1]     # y(x)  at x=x+h
    yp  = a21*yn[0] + a22*yn[1]     # y'(x) at x=x+h
    x  += h                         # x     at x=x+h
    return ((y,yp,x))

anz = 1000
h=2.*pi/anz
y0  =  1.    #inital value
y0p = -0.    #inital value 
x0  = 0.     #inital value 
loesung = []
y=(y0,y0p,x0)
for i in range(anz+1):
    yn = func(y,h)
    xn = yn[2]
    # print('{:8.4f} {:8.4f} {:8.4f}'.format(yn[0],yn[1],xn))
    y = yn
    loesung.append(y)

x=[v[2] for v in loesung]
y=[v[0] for v in loesung]
plot(x,y)
show()
