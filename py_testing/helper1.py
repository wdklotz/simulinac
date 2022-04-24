from math import sin,cos,degrees,radians
from matplotlib.pyplot import plot,show,legend,figure,subplot,axis
"""
    Funktion f(delta_phi) aus Anhang F in Trace3D manual
"""
def func(dphi):
    dp = radians(dphi)
    res = ((sin(dp)-dp*cos(dp))/dp)
    res = res*3./(dp*dp)
    res = res - sin(dp)/dp
    res = (15.*res)/(dp*dp)
    return res

print('Funktion 1-f(delta_phi) aus Anhang F in Trace3D manual')
x=[]
y=[]
dphi0 = 40.
for i in range(1,int(dphi0)):
    dphi = +(dphi0 - i)
    funcw = func(dphi)
    x.append(dphi)
    y.append(funcw)
    print('delta_phi[deg] {}\tf(delta_phi)[%] {:8.4f}'.format(dphi,(1.-funcw)*100.))
subplot(111).set_title('Funktion f(delta_phi) aus Anhang F in Trace3D manual')
pl=plot(x,y)
show()
    
