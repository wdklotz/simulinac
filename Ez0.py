import matplotlib.pyplot as plt
import numpy as np
from math import sin,cos,tan,pi,exp
from collections import namedtuple

from setutil import PARAMS,DEBUG,Proton

'''
Cavity E(z,r=0) field profile
'''
def NGauss(x,sig,mu):    # Gauss Normalverteilung
    res = exp(-(((x-mu)/sig)**2/2.))
    return res

def Kpoly(z,sigma,mu,E):
    """
    Calculate polynom coefficients
    """
    poly = []
    Interval = namedtuple('Interval',['zl','z0','zr','dz','b','a','E0','coeff'])
    anz = len(z)
    for i in range(0,anz-2,2):
        zl  = z[i]
        z0  = z[i+1]
        zr  = z[i+2]
        dz  =  z0-zl
        El  = E*NGauss(zl,sigma,mu)
        E0  = E*NGauss(z0,sigma,mu)
        Er  = E*NGauss(zr,sigma,mu)
        b = (Er+El-2*E0)/(2*E0*dz**2)
        a = (Er-El)/(2*E0*dz)
            # USE np.polyfit()
            # x   = np.array((zl,z0,zr))
            # y   = np.array((El,E0,Er))
            # coeff = np.polyfit(x,y,2)
            # b   = coeff[0]
            # a   = coeff[1]
            # E0  = coeff[2]
            # interval = Interval(zl,z0,zr,dz,b,a,E0,coeff)
            # ENDUSE np.polyfit()
        interval = Interval(zl,z0,zr,dz,b,a,E0,0.)
        poly.append(interval)
    return poly

def Epoly(z,poly):
    """
    eval polynomial fit
    """
    ix = -1
    for i in range(len(poly)):
        zl = poly[i].zl
        zr = poly[i].zr
        if zl <= z and z <= zr:
            ix = i
            break
    if ix <0:
        raise RuntimeError('arg out of range! {}'.format(z))
        # USE np.poly1d()
        # coeff = poly[ix].coeff
        # res   = np.poly1d(coeff)
        # res   = res(z)
        # ENDUSE np.poly1d()
    ival = poly[ix]
    z0 = ival.z0
    dz = z-z0
    E0 = ival.E0
    a  = ival.a
    b  = ival.b
    res = E0*(1.+a*dz+b*dz**2)
    return res

def V0(poly):
    v0 = []
    for i in range(len(poly)):
        E0 = poly[i].E0
        b  = poly[i].b
        dz = poly[i].dz
        v  = E0*(2*dz+2./3.*b*dz**3)
        v0.append(v)
    return v0

def Tk(poly,k):
    tk = []
    for i in range(len(poly)):
        a  = poly[i].a
        b  = poly[i].b
        dz = poly[i].dz
        rs = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
        rs = rs * (1.+b*dz**2-2*b/k**2+2*b*dz/k/tan(k*dz))
        tk.append(rs)
    return tk

def Sk(poly,k):
    sk = []
    for i in range(len(poly)):
        a  = poly[i].a
        b  = poly[i].b
        dz = poly[i].dz
        rs = 2*a*sin(k*dz)/k**2/(2*dz+2./3.*b*dz**3)
        rs = rs * (1.-k*dz/tan(k*dz))
        sk.append(rs)
    return sk

def Tkp(poly,k):
    tkp = []
    for i in range(len(poly)):
        a  = poly[i].a
        b  = poly[i].b
        dz = poly[i].dz
        rs = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
        rs = rs * ((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        tkp.append(rs)
    return tkp

def Skp(poly,k):
    skp = []
    for i in range(len(poly)):
        a  = poly[i].a
        b  = poly[i].b
        dz = poly[i].dz
        rs = 2*a*sin(k*dz)/k/(2*dz+2./3.*b*dz**2)
        rs = rs * (dz**2-2./k**2+2*dz/k/tan(k*dz))
        skp.append(rs)
    return skp

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

def display1(table,legend):
    z   = [+float(x[0]) for x in table]
    Ez  = [+float(x[2]) for x in table]
    plt.plot(z,Ez,label=legend)

def test0():
    '''
    Superfish data
    '''
    gap = 4.4
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
    display(Ez0_tab,'sf')
    # for x in Ez0_tab:
    #     print('z {}[cm]\tR {}[cm]\tEz(z,R) {}[MV/m]'.format(x[0],x[1],x[2]))

def test1():
    '''
    Gauss'che Normalverteilung
    '''
    gap = 4.4
    z = np.arange(0.,gap,0.044)
    sigm = 1.14
    Ez0_tab = [(x,0.,NGauss(x,sigm,0.)) for x in z]
    display(Ez0_tab,'NG')

    # Ez1_tab = [(x,0.,NGauss(x,sigm,gap)) for x in z]
    # Ez2_tab = []
    # for i in range(len(z)):
    #     x = z[i]
    #     Ez2 = Ez0_tab[i][2]+Ez1_tab[i][2]
    #     Ez2_tab.append((x,0.,Ez2))
    # cavlen = 2.5
    # print('sigma {}[cm], half-cavity length = {}[cm] = {}sigma'.
    #     format(sigm,cavlen,cavlen/sigm))
    # for x in Ez0_tab:
    #     print('z {}[cm]\tR {}[cm]\tEz(z,R) {}[MV/m]'.format(x[0],x[1],x[2]))

def test2():
    '''
    Second order polynomial fit
    '''
    particle = Proton(tkin=100.)
    beta     = particle.beta
    c        = PARAMS['lichtgeschwindigkeit']
    freq     = PARAMS['frequenz']
    k        = 2*pi*freq/(c*beta)

    anz   = 8            # nboff slices
    gap   = 4.4          #[cm] full gap length
    zl    = -2.*gap/2.   #left  interval boundary
    zr    = +2.*gap/2.   #right interval boundary
    sigma = gap/2./1.89  # sigma of NGauss (best fit with SF)
    # sigma = gap/2./2.2   # sigma of NGauss (best fit with SF)
    E0    = 1.           # top of NGauss   (best fit with SF)

    z = np.linspace(zl,zr,2*anz+1)
    # DEBUG('z',z)
    Ez0_tab = [(x,0.,E0*NGauss(x,sigma,0.)) for x in z]
    # display1(Ez0_tab,'slice')
    poly  = Kpoly(z,sigma,0.,E0*1.)
    # DEBUG('poly',poly)
    zstep = (zr-zl)/1000.
    z = np.arange(zl,zr,zstep)
    Ez0_tab = [(x,0.,Epoly(x, poly)) for x in z]
    display1(Ez0_tab,'poly')

    # poly1 = Kpoly(z,sigma,gap,E0*1.)  # next cavity
    # Ez1_tab = [(x,0.,Epoly(x,poly1)) for x in z]
    # display1(Ez1_tab,'poly1')
    # Ez2_tab = []
    # for i in range(len(z)):
    #     x = z[i]
    #     Ez2 = Ez0_tab[i][2]+Ez1_tab[i][2]
    #     Ez2_tab.append((x,0.,Ez2))
    # display1(Ez2_tab,'poly2')

    # TTF calculations
    v0  = V0(poly)
    tk  = Tk(poly,k)
    sk  = Sk(poly,k)
    tkp = Tkp(poly,k)
    skp = Skp(poly,k)
    DEBUG('V0',v0)
    DEBUG('T(k)',tk)
    DEBUG("T'(k)",tkp)
    DEBUG('S(k)',sk)
    DEBUG("S'(k)",skp)

if __name__ == '__main__':
    input_file = 'SF_WDK2g44.TBL'
    ax = pre_plt()
    test0()
    test1()
    test2()
    post_plt(ax)
        