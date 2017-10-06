import sys
sys.path.insert(0,'..')
from math import sqrt

from setutil import tblprnt,Proton

def fpsigma(psigma,soll):
    beta0 = soll.beta
    E0    = soll.e
    m0c2  = soll.e0
    res = (1+beta0**2*psigma)**2-(m0c2/E0)**2
    res = sqrt(res)/beta0-1.
    return res
def einsplusfpsigma(psigma,soll):
    return 1.+fpsigma(psigma,soll)
#conversion T3D ==> RipkenSchmidt (six)
def t3d2six(koord,soll,particle):
    x    = koord[0]
    xp   = koord[1]
    y    = koord[2]
    yp   = koord[3]
    z    = koord[4]
    dp2p = koord[5]

    E0    = soll.e
    beta0 = soll.beta
    m0c2  = soll.e0
    gb    = particle.gamma_beta
    beta  = particle.beta

    px    = gb*m0c2/E0*xp
    py    = gb*m0c2/E0*yp
    psigma = ((beta0/beta/(1.-dp2p))-1.)/beta0**2
    return (x,px,y,py,z,psigma)
# conversion RipkenSchmidt (six) ==> T3D
def six2t3d(koord,soll,particle):
    x      = koord[0]
    px     = koord[1]
    y      = koord[2]
    py     = koord[3]
    sigma  = koord[4]
    psigma = koord[5]

    beta  = particle.beta
    gb    = particle.gamma_beta
    beta0 = soll.beta
    E0    = soll.e
    m0c2  = soll.e0

    xp   = px/(gb*m0c2/E0)
    yp   = py/(gb*m0c2/E0)
    z    = sigma
    dp2p = 1.-beta0/beta/(1.+beta0**2*psigma)
    return (x,xp,y,yp,z,dp2p)
# Ripken-Schnidt (six) map
def sixmap(koord,soll,particle,l):
    xi       = koord[0]
    pxi      = koord[1]
    yi       = koord[2]
    pyi      = koord[3]
    sigmai   = koord[4]
    psigmai  = koord[5]

    beta0    = soll.beta
    beta     = particle.beta
    xf       = xi + pxi/einsplusfpsigma(psigmai,soll)*l
    pxf      = pxi
    yf       = yi + pyi/einsplusfpsigma(psigmai,soll)*l
    pyf      = pyi
    sigmaf  = sigmai + (1.-(beta0/beta)*(1.+0.5*(pxi**2+pyi**2)/einsplusfpsigma(psigmai,soll)**2))*l
    psigmaf = psigmai
    return(xf,pxf,yf,pyf,sigmaf,psigmaf)
# T3D matrix map
def t3dmap(koord,particle,l):
    xi    = koord[0]
    xpi   = koord[1]
    yi    = koord[2]
    ypi   = koord[3]
    zi    = koord[4]
    dp2pi = koord[5]

    gamma = particle.gamma

    xf     = xi+l*xpi
    xpf    = xpi
    yf     = yi+l*ypi
    ypf    = ypi
    zf     = zi+l/gamma**2*dp2pi
    dp2pf  = dp2pi
    return(xf,xpf,yf,xpf,zf,dp2pf)
##body 
header_six=[
        'x',
        "px",
        'y',
        "py",
        'sigma',
        'psigma'
        ]
header_t3d=[
        'x',
        "x'",
        'y',
        "y'",
        'z',
        'dp/p'
        ]
tkin0      = 100.
soll       = Proton(tkin=tkin0)
E0         = soll.e          # E-soll [MeV]
p0         = soll.p          # cp-soll [MeV]
m0c2       = soll.e0         # rest mass soll [Mev]

dp2pi    = 1.e-2
p        = p0/(1.-dp2pi)
E        = sqrt(p**2+m0c2**2) #E aus dp2p und p0
tkini    = E-m0c2
particle = Proton(tkin=tkini)

l          = 50.e-3          #[m]

# x,xp,y,yp,z,dp/p  T3D Kordinaten
xi  = yi  = 1.e-3
xpi = ypi = 1.e-3
zi  = 1.e-3
# dp2pi = (p-p0)/p
koord0 = (xi,xpi,yi,ypi,zi,dp2pi)
row = [[
        '{:9.6e}'.format(xi),
        '{:9.6e}'.format(xpi),
        '{:9.6e}'.format(yi),
        '{:9.6e}'.format(ypi),
        '{:9.6e}'.format(zi),
        '{:9.6e}'.format(dp2pi)
        ]]
print('T3D IN\n'+tblprnt(header_t3d,row))

# x,px,y,py,sg,psg   kanonische Koordinaten
koord_six = t3d2six(koord0,soll,particle)
xi      = koord_six[0]
pxi     = koord_six[1] #px
yi      = koord_six[2]
pyi     = koord_six[3] #py
sigmai  = koord_six[4] #sigma
psigmai = koord_six[5] #psigma
row = [[
        '{:10.6g}'.format(xi),
        '{:10.6g}'.format(pxi),
        '{:10.6g}'.format(yi),
        '{:10.6g}'.format(pyi),
        '{:10.6g}'.format(sigmai),
        '{:10.6g}'.format(psigmai)
        ]]
print('six IN \n'+tblprnt(header_six,row))

koord = six2t3d(koord_six,soll,particle)
xr    = koord[0]
xpr   = koord[1]
yr    = koord[2]
ypr   = koord[3]
zr    = koord[4]
dp2pr = koord[5]
row = [[
        '{:9.6e}'.format(xr),
        '{:9.6e}'.format(xpr),
        '{:9.6e}'.format(yr),
        '{:9.6e}'.format(ypr),
        '{:9.6e}'.format(zr),
        '{:9.6e}'.format(dp2pr)
        ]]
print('six IN zurueck zu t3d IN\n'+tblprnt(header_t3d,row))

koord_six_out = sixmap(koord_six,soll,particle,l)
xf       = koord_six_out[0]
pxf      = koord_six_out[1]
yf       = koord_six_out[2]
pyf      = koord_six_out[3]
sigmaf   = koord_six_out[4]
psigmaf  = koord_six_out[5]
row = [[
        '{:10.6g}'.format(xf),
        '{:10.6g}'.format(pxf),
        '{:10.6g}'.format(yf),
        '{:10.6g}'.format(pyf),
        '{:10.6g}'.format(sigmaf),
        '{:10.6g}'.format(psigmaf)
        ]]
print('six OUT \n'+tblprnt(header_six,row))

koord = six2t3d(koord_six_out,soll,particle)
xf    = koord[0]
xpf   = koord[1]
yf    = koord[2]
ypf   = koord[3]
zf    = koord[4]
dp2pf = koord[5]
row = [[
        '{:9.6e}'.format(xf),
        '{:9.6e}'.format(xpf),
        '{:9.6e}'.format(yf),
        '{:9.6e}'.format(ypf),
        '{:9.6e}'.format(zf),
        '{:9.6e}'.format(dp2pf)
        ]]
print('six OUT in T3D Koordinaten\n'+tblprnt(header_t3d,row))

koord_t3d = t3dmap(koord0,particle,l)
xf     = koord_t3d[0]
xpf    = koord_t3d[1]
yf     = koord_t3d[2]
ypf    = koord_t3d[3]
zf     = koord_t3d[4]
dp2pf  = koord_t3d[5]
row = [[
        '{:9.6e}'.format(xf),
        '{:9.6e}'.format(xpf),
        '{:9.6e}'.format(yf),
        '{:9.6e}'.format(ypf),
        '{:9.6e}'.format(zf),
        '{:9.6e}'.format(dp2pf)
        ]]
print('T3D OUT\n'+tblprnt(header_t3d,row))


