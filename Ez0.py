import matplotlib.pyplot as plt
import numpy as np
from math import sin,cos,tan,pi,exp

from setutil import PARAMS,DEBUG,Proton

'''
Cavity E(z,r=0) field profile
'''
def NGauss(x,sig,mu):    # Gauss Normalverteilung
    res = exp(-(((x-mu)/sig)**2/2.))
    return res

def Epoly(z,poly):
    """
    Plynomial fit
    """
    ix = -1
    for i in range(len(poly)):
        zl = poly[i][0]
        z0 = poly[i][1]
        zr = poly[i][2]
        if zl <= z and z <= zr:
            ix = i
            break
    if ix <0:
        raise RuntimeError('arg out of range! {}'.format(z))
    E0 = poly[ix][4]
    a  = poly[ix][5]
    b  = poly[ix][6]
    dz = z-z0
    res = E0*(1.+ a*dz+b*dz**2)
    return res

def V0(poly):
    v0 = []
    for i in range(len(poly)):
        zl = poly[i][0]
        z0 = poly[i][1]
        E0 = poly[i][4]
        b  = poly[i][6]
        dz = z0-zl
        v  = E0*(2*dz+2./3.*b*dz**3)
        v0.append(v)
    return v0

def Tk(poly,k):
    tk = []
    for i in range(len(poly)):
        zl = poly[i][0]
        z0 = poly[i][1]
        a  = poly[i][5]
        b  = poly[i][6]
        dz = z0-zl
        rs = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
        rs = rs * (1.+b*dz**2-2*b/k**2+2*b*dz/k/tan(k*dz))
        tk.append(rs)
    return tk

def Sk(poly,k):
    sk = []
    for i in range(len(poly)):
        zl = poly[i][0]
        z0 = poly[i][1]
        a  = poly[i][5]
        b  = poly[i][6]
        dz = z0-zl
        rs = 2*a*sin(k*dz)/k**2/(2*dz+2./3.*b*dz**3)
        rs = rs * (1.-k*dz/tan(k*dz))
        sk.append(rs)
    return sk

def Tkp(poly,k):
    tkp = []
    for i in range(len(poly)):
        zl = poly[i][0]
        z0 = poly[i][1]
        a  = poly[i][5]
        b  = poly[i][6]
        dz = z0-zl
        rs = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
        rs = rs * ((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        tkp.append(rs)
    return tkp

def Skp(poly,k):
    skp = []
    for i in range(len(poly)):
        zl = poly[i][0]
        z0 = poly[i][1]
        a  = poly[i][5]
        b  = poly[i][6]
        dz = z0-zl
        rs = 2*a*sin(k*dz)/k/(2*dz+2./3.*b*dz**2)
        rs = rs * (dz**2-2./k**2+2*dz/k/tan(k*dz))
        skp.append(rs)
    return skp

def Kpoly(z,sigma,E):
    """
    Calculate polynom coefficients
    """
    interval = []
    anz = int((len(z)-1)/2)
    for i in range(0,anz+2,2):
        zl  = z[i]
        z0  = z[i+1]
        zr  = z[i+2]
        dz  =  z0-zl
        El  = E*NGauss(zl,sigma,0.)
        E0  = E*NGauss(z0,sigma,0.)
        Er  = E*NGauss(zr,sigma,0.)
        a   = (Er-El)/(2*E0*dz)
        b   = (Er+El-2*E0)/(2*E0*dz**2)
        interval.append((zl,z0,zr,dz,E0,a,b))
    return interval

def pre_plt():
    ax  = plt.subplot(111)
    ax.set_title(input_file)
    ax.set_ylabel('Ez0 [MV/m]')
    ax.set_xlabel('z [cm]')
    return ax

def post_plt(ax):
    plt.legend(loc='lower right',fontsize='x-small')
    plt.show()

def display(table,legend):
    zp   = [+float(x[0]) for x in table]
    Ezp  = [+float(x[2]) for x in table]
    zn   = [-float(x[0]) for x in reversed(table)]
    Ezn  = [+float(x[2]) for x in reversed(table)]
    plt.plot(zn+zp,Ezn+Ezp,label=legend)

def test0():
    '''
    Superfish data
    '''
    Ez0_tab = []
    Ez_max = 1.
    with open(input_file,'r') as f:
        lines = list(f)
        for cnt,line in enumerate(lines[41:-2]):
            # print(line,end='')
            stripped    = line.strip()
            (z,sep,aft) = stripped.partition(' ')
            stripped    = aft.strip()
            (R,sep,aft) = stripped.partition(' ')
            stripped     = aft.strip()
            (Ez,sep,aft) = stripped.partition(' ')
            if cnt == 0:
                Ez_max = float(Ez)      # normalization factor
            Ez = float(Ez)/Ez_max
            Ez0_tab.append((z,R,Ez))
    # for x in Ez0_tab:
    #     print('z {}[cm]\tR {}[cm]\tEz(z,R) {}[MV/m]'.format(x[0],x[1],x[2]))
    display(Ez0_tab,'sf')

def test1():
    '''
    Gauss'che Normalverteilung
    '''
    z = np.arange(0.,4.4,0.044)
    sigm = 1.14
    Ez0_tab = [(x,0.,NGauss(x,sigm,0.)) for x in z]
    cavlen = 2.5
    # print('sigma {}[cm], half-cavity length = {}[cm] = {}sigma'.
    #     format(sigm,cavlen,cavlen/sigm))
    # for x in Ez0_tab:
    #     print('z {}[cm]\tR {}[cm]\tEz(z,R) {}[MV/m]'.format(x[0],x[1],x[2]))
    display(Ez0_tab,'NG')

def test2():
    '''
    Second order polynomial fit
    '''
    particle = Proton(tkin=100.)
    beta     = particle.beta
    c        = PARAMS['lichtgeschwindigkeit']
    freq     = PARAMS['frequenz']
    k        = 2*pi*freq/(c*beta)

    anz   = 3            # nboff slices
    gap   = 4.4          #[cm] full gap length
    zl    = 0.           #left  interval boundary
    zr    = gap/2.       #right interval boundary
    sigma = (zr-zl)/1.89 # sigma of NGauss (best fit with SF)
    E0    = 1.0          # top of NGauss   (best fit with SF)

    zstep = (zr-zl)/(2*anz)
    z = np.arange(zl,zr+zstep,zstep)
    # DEBUG('z',z)
    # Ez0_tab = [(x,0.,E0*NGauss(x,sigma,0.)) for x in z]
    # display(Ez0_tab,'slice')
    poly = Kpoly(z,sigma,E0)
    # DEBUG('poly',poly)

    zstep = (zr-zl)/1000.
    z = np.arange(zl,zr,zstep)
    Ez0_tab = [(x,0.,Epoly(x,poly)) for x in z]
    display(Ez0_tab,'poly')

    # TTF calculations
    v0 = V0(poly)
    DEBUG('v0',v0)
    tk = Tk(poly,k)
    DEBUG('T(k)',tk)
    sk = Sk(poly,k)
    DEBUG('S(k)',sk)
    tkp = Tkp(poly,k)
    DEBUG("T'(k)",tkp)
    skp = Skp(poly,k)
    DEBUG("S'(k)",skp)

if __name__ == '__main__':
    input_file = 'SF_WDK2g44.TBL'
    ax = pre_plt()
    test0()
    test1()
    test2()
    post_plt(ax)
        