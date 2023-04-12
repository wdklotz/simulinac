import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from math import sqrt

a     = 2
b     = 1./2.
eps   = a * b
alfa  = -15./8. 
beta  = 17./8.
gamma = 17./8.
vLines = [0.686,a/sqrt(2),sqrt((a**2+b**2)/2)]
elli0 = Ellipse((0,0),2*a,2*a,fill=False,color='green')
elli1 = Ellipse((0,0),2*b,2*a,fill=False,color='blue')
elli2 = Ellipse((0,0),2*b,2*a,angle=-45,fill=False,color='red')
# elli3 = Ellipse((0,0),2*0.5,2*2.0,angle=-45,fill=False,color='red')

xdata = np.linspace(-1.6,1.6, num=80)
# print(xdata)
def make_elli(alfa,beta,gamma,eps,x):
    pts=[]
    for x in xdata:
        crit = (alfa*x)**2 - beta*(gamma*x**2-eps)
        if crit >= 0.:
            sqcrit = sqrt(4.*crit)
            yp = (-2*alfa*x + sqcrit)/(2*beta)
            ym = (-2*alfa*x - sqcrit)/(2*beta)
            pp = (x,yp)
            pm = (x,ym)
            pts.append(pp)
            pts.append(pm) 
    return pts

fig,ax = plt.subplots(figsize=(8,8))
ax.add_patch(elli0)
ax.add_patch(elli1)
ax.add_patch(elli2)

pts = make_elli(alfa,gamma,beta,eps,xdata)
xc=[x[0] for x in pts]
yc=[y[1] for y in pts]
ax.plot(xc,yc,'x',markersize=4,color='black')

pts = make_elli(-3.5/2.,2.01,2.01,eps,xdata)
xc=[x[0] for x in pts]
yc=[y[1] for y in pts]
ax.plot(xc,yc,'x',markersize=4,color='sienna')

ax.vlines(vLines,-2.,+2.)
plt.xlim(-2.,+2.) 
plt.ylim(-2.,+2.)
plt.grid()
plt.show()

