import sys
sys.path.insert(0,'..')

from math import sqrt
from setutil import tblprnt,Proton,PARAMS

def f(psigma,beta0,E0):
    res = (1+beta0**2*psigma)**2-(m0c2/E0)**2
    res = sqrt(res)/beta0-1.
    return res
#conversion T3D ==> RipkenSchmidt
def t3d2rip(koord,soll,particle):
    x    = koord[0]
    xp   = koord[1]
    y    = koord[2]
    yp   = koord[3]
    z    = koord[4]
    dp2p = koord[5]
    E0    = soll.e
    beta0 = soll.beta
    gb    = particle.gamma_beta
    beta  = particle.beta
    gamma = particle.gamma
    px    = gb*m0c2/E0*xp
    py    = gb*m0c2/E0*yp
    psigma = (beta/beta0)**2*gamma*dp2p
    return (x,px,y,py,z,psigma)
# conversion RipkenSchmidt ==> T3D
def rip2t3d(koord,soll,particle):
    beta  = particle.beta
    gb    = particle.gamma_beta
    gamma = particle.gamma
    beta0 = soll.beta
    E0    = soll.e
    xf      = koord[0]
    pxf     = koord[1]
    yf      = koord[2]
    pyf     = koord[3]
    sigmaf  = koord[4]
    psigmaf = koord[5]
    xpf   = pxf/(gb*m0c2/E0)
    ypf   = pyf/(gb*m0c2/E0)
    zf    = sigmaf
    dp2pf = psigmaf/(((beta/beta0)**2*gamma))
    return (xf,xpf,yf,ypf,zf,dp2pf)
# Ripken-Schnidt map
def ripmap(koord,soll,particle,l):
    beta0    = soll.beta
    E0       = soll.e
    beta     = particle.beta
    xi       = koord[0]
    pxi      = koord[1]
    yi       = koord[2]
    pyi      = koord[3]
    sigmai   = koord[4]
    psigmai  = koord[5]
    xf   = xi + pxi/(1+f(psigmai,beta0,E0))*l
    pxf  = pxi
    yf   = yi + pyi/(1+f(psigmai,beta0,E0))*l
    pyf  = pyi
    sigmaf  = sigmai + (1.-(beta0/beta)*(1.+0.5*(pxi**2+pyi**2)/(1+f(psigmai,beta0,E0))**2))*l
    psigmaf = psigmai
    return(xf,pxf,yf,pyf,sigmaf,psigmaf)
# T3D matrix map
def t3dmap(koord,particle,l):
    gamma = particle.gamma
    xi    = koord[0]
    xpi   = koord[1]
    yi    = koord[2]
    ypi   = koord[3]
    zi    = koord[4]
    dp2pi = koord[5]
    xf     = xi+l*xpi
    xpf    = xpi
    yf     = yi+l*ypi
    ypf    = ypi
    zf     = zi+l/gamma**2*dp2pi
    dp2pf  = dp2pi
    return(xf,xpf,yf,xpf,zf,dp2pf)
    
m0c2 = PARAMS['proton_mass']
tkin0      = 100.
soll       = Proton(tkin=tkin0)
E0         = soll.e          # E-soll

eta        = 1.e-2           #[%] dE/E0=(E-E0)/E0
E          = (eta+1)*E0
tkini      = E-m0c2
particle   = Proton(tkin=tkini)
l          = 50.e-3          #[m]

# x,xp,y,yp,z,dp/p  T3D Kordinaten
xi  = yi  = 1.e-3
xpi = ypi = 1.e-3
zi  = dp2pi = 1.e-3
header=[
        'x',
        "x'",
        'y',
        "y'",
        'z',
        'dp/p'
        ]
row = [[
        '{:9.6e}'.format(xi),
        '{:9.6e}'.format(xpi),
        '{:9.6e}'.format(yi),
        '{:9.6e}'.format(ypi),
        '{:9.6e}'.format(zi),
        '{:9.6e}'.format(dp2pi)
        ]]
print('T3D IN\n'+tblprnt(header,row))

# x,px,y,py,sg,psg   kanonische Koordinaten
koord = t3d2rip((xi,xpi,yi,ypi,zi,dp2pi),soll,particle)
xi      = koord[0]
pxi     = koord[1] #px
yi      = koord[2]
pyi     = koord[3] #py
sigmai  = koord[4] #sigma
psigmai = koord[5] #psigma
header=[
        'x',
        "px",
        'y',
        "py",
        'sigma',
        'psigma'
        ]
row = [[
        '{:10.6g}'.format(xi),
        '{:10.6g}'.format(pxi),
        '{:10.6g}'.format(yi),
        '{:10.6g}'.format(pyi),
        '{:10.6g}'.format(sigmai),
        '{:10.6g}'.format(psigmai)
        ]]
print('RIP IN \n'+tblprnt(header,row))

koord = ripmap((xi,pxi,yi,pyi,sigmai,psigmai),soll,particle,l)
xf       = koord[0]
pxf      = koord[1]
yf       = koord[2]
pyf      = koord[3]
sigmaf   = koord[4]
psigmaf  = koord[5]
row = [[
        '{:10.6g}'.format(xf),
        '{:10.6g}'.format(pxf),
        '{:10.6g}'.format(yf),
        '{:10.6g}'.format(pyf),
        '{:10.6g}'.format(sigmaf),
        '{:10.6g}'.format(psigmaf)
        ]]
print('RIP OUT \n'+tblprnt(header,row))

koord = rip2t3d((xf,pxf,yf,pyf,sigmaf,psigmaf),soll,particle)
xf    = koord[0]
xpf   = koord[1]
yf    = koord[2]
ypf   = koord[3]
zf    = koord[4]
dp2pf = koord[5]
header=[
        'x',
        "x'",
        'y',
        "y'",
        'z',
        'dp/p'
        ]
row = [[
        '{:9.6e}'.format(xf),
        '{:9.6e}'.format(xpf),
        '{:9.6e}'.format(yf),
        '{:9.6e}'.format(ypf),
        '{:9.6e}'.format(zf),
        '{:9.6e}'.format(dp2pf)
        ]]
print('RIP OUT in T3D Koordinaten\n'+tblprnt(header,row))

koord = t3dmap((xi,xpi,yi,ypi,zi,dp2pi),particle,l)
xf     = koord[0]
xpf    = koord[1]
yf     = koord[2]
ypf    = koord[3]
zf     = koord[4]
dp2pf  = koord[5]
row = [[
        '{:9.6e}'.format(xf),
        '{:9.6e}'.format(xpf),
        '{:9.6e}'.format(yf),
        '{:9.6e}'.format(ypf),
        '{:9.6e}'.format(zf),
        '{:9.6e}'.format(dp2pf)
        ]]
print('T3D OUT\n'+tblprnt(header,row))


