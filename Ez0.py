import matplotlib.pyplot as plt
import numpy as np
from math import sin,cos,tan,pi,exp,fabs,pow,sqrt,fmod
from collections import namedtuple

from setutil import PARAMS,DEBUG,Proton

'''
Cavity E(z,r=0) field profile
'''
def NGauss(x,sig,mu):    # Gauss Normalverteilung
    res = exp(-(((x-mu)/sig)**2/2.))
    return res

def SFdata(input_file,Epeak=1.):
    '''
    Superfish data
    '''
    Dpoint = namedtuple('DataPoint',['z','R','Ez'])  # data-point structure
    zp = []
    rp = []
    ep = []
    with open(input_file,'r') as f:
        lines = list(f)
        for line in lines[41:-2]:
            # print(line,end='')
            stripped    = line.strip()
            (z,sep,aft) = stripped.partition(' ')
            z =float(z)
            stripped    = aft.strip()
            (R,sep,aft) = stripped.partition(' ')
            R = float(R)
            stripped     = aft.strip()
            (Ez,sep,aft) = stripped.partition(' ')
            Ez = float(Ez)
            zp.append(z)
            rp.append(R)
            ep.append(Ez)
    zprev = [-x for x in reversed(zp[1:])]
    zp = zprev+zp
    rprev = [-x for x in reversed(rp[1:])]
    rp = rprev+rp
    eprev = [x for x in reversed(ep[1:])]
    ep =eprev+ep
    emax = max(ep)
    enorm = Epeak/emax
    Ez0_tab = []
    for i in range(len(zp)):   # normalize and pack
        ep[i] = ep[i]*enorm
        dpoint = Dpoint(zp[i],rp[i],ep[i])
        Ez0_tab.append(dpoint)
    return Ez0_tab

            # EXPERIMENTAL
# b0 = 2.5052367      from Abramowitz-Stegun
# b2 = 1.2831204
# b4 = 0.2264718
# b6 = 0.1306469
# b8 = -0.0202490
# b10= 0.0039132
# zcof = [b0,0,b2,0,b4,0,b6,0,b8,0,b10]

# def Z(x,sigma,mu):       #from Abramowitz-Stegun
#     x =(x-mu)/sigma
#     b0 = zcof[0]
#     b2 = zcof[2]
#     b4 = zcof[4]
#     b6 = zcof[6]
#     b8 = zcof[8]
#     b10= zcof[10]
#     res = sqrt(2*pi)/(b0+b2*pow(x,2)+b4*pow(x,4)+b6*pow(x,6)+b8*pow(x,8)+b10*pow(x,10))
#     return res

# def ZR():                  # Verusch mit numpy library
#     quo,rem = np.polynomial.polynomial.polydiv([1.,0.,0.,0.,0.,0.],zcof)
#     print('quo,rem',quo,rem)

# def Zinv(x,sigma,mu):       # Versuch mit reziprokem Polynom
#     b = zcof
#     c = [0. for x in range(11)]
#     c[0] = 1./b[0]
#     for n in range(1,11,1):
#         for i in range(0,n,1):
#             print('c[{}],b[{}] {} {}'.format(i,n-i,c[i],b[n-i]))
#             c[n] -= (c[i]*b[n-i])/b[0]
#     print(c)
#     x = (x-mu)/sigma
#     b0 = c[0]
#     b2 = c[2]
#     b4 = c[4]
#     b6 = c[6]
#     b8 = c[8]
#     b10= c[10]
#     res = sqrt(2*pi)*(b0+b2*pow(x,2)+b4*pow(x,4)+b6*pow(x,6)+b8*pow(x,8)+b10*pow(x,10))
#     return res

# def polint(xa,ya,n,x): # Versuch Uebersetzung von Routinen von Numerical recipies
#     """
#     Given arrays xa,ya and given a value x, this routine returns a value y,
#     and an error estimate dy. If P(x) is the polynomial of degre n-1 such that
#     P(xai) = yai, i = (1,...,n), then the returned value y = P(x).
#     Numerical recipies in C, Cambridge Univ. Press, 1988
#     """
#     def IX(n):
#         return n-1
# 
#     # DEBUG('xa',xa)
#     # DEBUG('ya',ya)
#     # DEBUG('n ',n)
#     ns = 1
#     dif = fabs(x-xa[IX(1)])
#     c = [ya[i] for i in range(n)]
#     d = [ya[i] for i in range(n)]
#     for i in range(1,n+1):
#         dift = fabs(x-xa[IX(i)])
#         # DEBUG('dift ',dift)
#         if dift < dif:
#             ns = i
#             dif = dift
#         c[IX(i)] = ya[IX(i)]
#         d[IX(i)] = ya[IX(i)]
#     # DEBUG('ns ',ns)
#     ns -= 1
#     y = ya[IX(ns)]
#     # DEBUG('ns ',ns)
#     # DEBUG('y ',y)
#     dy = 0.
#     for m in range(1,n):
#         # DEBUG('range(1,n-m+1) ',range(1,n-m+1))
#         for i in range(1,n-m+1):
#             ho = xa[IX(i)]-x
#             hp = xa[IX(i+m)]-x
#             w = c[IX(i+1)]-d[IX(i)]
#             den = ho-hp
#             # DEBUG('m,i,ho,hp,w,den {},{},{},{},{},{}'.format(m,i,ho,hp,w,den))
#             if den == 0.:
#                 raise RuntimeError('Error in routine POLINT')
#             den = w/den
#             d[IX(i)] = hp*den
#             c[IX(i)] = ho*den
#         if (2*ns) < (n-m):
#             dy = c[IX(ns+1)]
#         else:
#             ns -= 1
#             dy  = d[IX(ns)]
#         y += dy
#     return y,dy

# def polcof(x,y,n):
#     """
#     Given arrays x,y containing a tabulated function f=f(x), this routine
#     returns an array of coefficients cof[i], such that y(xi)=SUMi(cof[i]*xi**i)
#     Numerical recepies in C, Cambridge Univ. Press, 1988
#     """
#     xw = [x[i] for i in range(len(x))]
#     yw = [y[i] for i in range(len(y))]
#     cof = np.zeros(len(x))
# 
#     for j in range(n+1):
#         # DEBUG('j ',j)
#         cof[j],dy = polint(xw,yw,n+1-j,0.0)
#         # DEBUG('cof[j]',cof[j])
#         xmin = 1.0e38
#         k = -1
#         for i in range(n-j+1):
#             if fabs(xw[i]) < xmin:
#                 xmin = fabs(xw[i])
#                 k = i
#             if xw[i] != 0.:
#                 yw[i] = (yw[i]-cof[j])/xw[i]
#         for i in range(k+1,n-j+1):
#             yw[i-1] = yw[i]
#             xw[i-1] = xw[i]
#     return cof
            # END-EXPERIMENTAL

def KpolySF(Ez0_tab,anz=4):
    """
    Calculate polynom coefficients from SF data
    """
    def indexer(nbslices,M):
        """
            nbslices = nboff slices
            N = nboff half-length-slices
            M = nboff poly-points
            n = nboff poly-intervals/half-length-slice
        """
        N=2*nbslices    # factor 2 more intervals than slices
        if N>M: 
            raise RuntimeError('nboff slices must be <= {}'.format(int((M-1)/2)))
        M = int(M-fmod(M,N))
        n = int(M/N)
        # DEBUG('(N,M,n)=({},{},{})'.format(N,M,n))
        for i in range(0,M,2*n):
            # DEBUG('(i,i+n,i+2*n)={},{},{}'.format(i,i+n,i+2*n))
            yield((i,i+n,i+2*n))
    
    poly = []
    Interval = namedtuple('Interval',['zl','z0','zr','dz','b','a','E0','coeff'])
    for (il,i0,ir) in indexer(anz,len(Ez0_tab)):
        # DEBUG('(il,i0,ir) ',((il,i0,ir)))
        zl = Ez0_tab[il].z
        z0 = Ez0_tab[i0].z
        zr = Ez0_tab[ir].z
        El = Ez0_tab[il].Ez
        E0 = Ez0_tab[i0].Ez
        Er = Ez0_tab[ir].Ez
        dz = z0-zl
        b  = (Er+El-2*E0)/(2*E0*dz**2)   # Langrange 3 Punkt Interpolation 
        a  = (Er-El)/(2*E0*dz)           # getestet mit Bleistift u. Papier
        interval = Interval(zl,z0,zr,dz,b,a,E0,0.)
        poly.append(interval)
    return poly

def Kpoly(z,sigma,mu,E):
    """
    Calculate polynom coefficients from NG formula
    """
    poly = []
    Interval = namedtuple('Interval',['zl','z0','zr','dz','b','a','E0','coeff'])
    anz = len(z)
    for i in range(0,anz-2,2):
        zl  = z[i]
        z0  = z[i+1]
        zr  = z[i+2]
        dz  = z0-zl
        El  = E*NGauss(zl,sigma,mu)
        E0  = E*NGauss(z0,sigma,mu)
        Er  = E*NGauss(zr,sigma,mu)
        b = (Er+El-2*E0)/(2*E0*dz**2)   # Langrange 3 Punkt Interpolation 
        a = (Er-El)/(2*E0*dz)           # getetstet mit Bleistift u. Papier
        interval = Interval(zl,z0,zr,dz,b,a,E0,0.)
        poly.append(interval)

            # USE np.polyfit() - liefert gleiches Resultat wie Langrange 3 Punkt
        # x   = np.array((zl,z0,zr))
        # y   = np.array((El,E0,Er))
        # coeff = np.polyfit(x,y,2)
        # b   = coeff[0]
        # a   = coeff[1]
        # E0  = coeff[2]
        # interval = Interval(zl,z0,zr,dz,b,a,E0,coeff)

        # x = [zl,z0,zr]
        # y = [El,E0,Er]
        # cof = polcof(x,y,2)
        # print('cof =====',cof)
        # print('E0,a,b ',E0,a,b,'\n')
            # END-USE np.polyfit()

    return poly

def Ipoly(z,poly):
    """
    interpolate with polynomial fit
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

            # USE np.polyfit()
            # liefert gleiches Resultat wie Langrange 3 Punkt
    # coeff = poly[ix].coeff
    # res   = np.poly1d(coeff)
    # res   = res(z)
            # END-USE np.poly1d()

    ival = poly[ix]
    z0 = ival.z0
    dz = z-z0
    E0 = ival.E0
    a  = ival.a
    b  = ival.b
    res = E0*(1.+a*dz+b*dz**2)
    return res

def V0(poly,zl,zr):
    v0 = []
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz  = zir-zil
        if zil < zl or zir > zr: continue
        E0 = poly[i].E0
        b  = poly[i].b
        dz = poly[i].dz
        v  = E0*(2*dz+2./3.*b*dz**3)
        v0.append(v)
    return v0

def Tk(poly,k,zl,zr):
    tk = []
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz  = zir-zil
        if zil < zl or zir > zr: continue
        # DEBUG('Tk: (i,dz,zl,zil,zir,zr)=({:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(i,dz,zl,zil,zir,zr))
        a  = poly[i].a
        b  = poly[i].b
        dz = poly[i].dz
        f1 = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
        f2 = 1.+b*dz**2-2.*b/k**2*(1.-k*dz/tan(k*dz))
        # DEBUG('Tk: (a,b,dz,f1,f2)={:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(a,b,dz,f1,f2))
        rs = f1*f2
        tk.append(rs)
    return tk

def Sk(poly,k,zl,zr):
    sk = []
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz  = zir-zil
        if zil < zl or zir > zr: continue
        a  = poly[i].a
        b  = poly[i].b
        dz = poly[i].dz
        rs = 2*a*sin(k*dz)/k**2/(2*dz+2./3.*b*dz**3)
        rs = rs * (1.-k*dz/tan(k*dz))
        sk.append(rs)
    return sk

def Tkp(poly,k,zl,zr):
    tkp = []
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz  = zir-zil
        if zil < zl or zir > zr: continue
        a  = poly[i].a
        b  = poly[i].b
        dz = poly[i].dz
        rs = 2*sin(k*dz)/k/(2*dz+2./3.*b*dz**3)
        rs = rs * ((1.+3*b*dz**2-6*b/k**2)/k-dz/tan(k*dz)*(1.+b*dz**2-6*b/k**2))
        tkp.append(rs)
    return tkp

def Skp(poly,k,zl,zr):
    skp = []
    for i in range(len(poly)):
        zil = poly[i].zl
        zir = poly[i].zr
        dz  = zir-zil
        if zil < zl or zir > zr: continue
        a  = poly[i].a
        b  = poly[i].b
        dz = poly[i].dz
        rs = 2*a*sin(k*dz)/k/(2*dz+2./3.*b*dz**2)
        rs = rs * (dz**2-2./k**2+2*dz/k/tan(k*dz))
        skp.append(rs)
    return skp

def pre_plt(input_file):
    ax  = plt.subplot(111)
    ax.set_ylabel('Ez0 [MV/m]')
    ax.set_xlabel('z [cm]')
    ax.set_title(input_file)
    return ax

def post_plt(ax):
    plt.legend(loc='lower right',fontsize='x-small')
    plt.show()

def displayLR(table,legend):
    zp   = [+float(x[0]) for x in table]
    Ezp  = [+float(x[2]) for x in table]
    zn   = [-float(x[0]) for x in reversed(table)]
    Ezn  = [+float(x[2]) for x in reversed(table)]
    plt.plot(zn+zp,Ezn+Ezp,label=legend)

def display(table,legend):
    z   = [+float(x[0]) for x in table]
    Ez  = [+float(x[2]) for x in table]
    plt.plot(z,Ez,label=legend)

def test0(input_file,ax,plot=True):
    Ez0_tab = SFdata(input_file)
    if plot: 
        display(Ez0_tab,'SF')
    
def test1():
    '''
    Gauss'che Normalverteilung
    '''
    gap = 4.4
    z = np.arange(0.,gap,gap/500.)
    sigma = 1.14
    Ez0_tab = [(x,0.,NGauss(x,sigma,0.)) for x in z]
    displayLR(Ez0_tab,'NG')

            # EXPERIMENTAL
    # EzZ_tab = [(x,0.,Z(x,sigma,0.01 )) for x in z]
    # display(EzZ_tab,'Z')
            # END-EXPERIMENTAL

            # OVERAP?
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
            # END-OVERAP?

def test2():
    '''
    Second order polynomial fit with Lagrange 3 point formula to NG formula
    '''
    particle = Proton(tkin=100.)
    beta     = particle.beta
    c        = PARAMS['lichtgeschwindigkeit']
    freq     = PARAMS['frequenz']
    k        = 2*pi*freq/(c*beta)

    anz   = 6           # nboff slices
    gap   = 4.8          #[cm] full gap length
    zl    = -gap/2.      #left  interval boundary
    zr    = gap/2.       #right interval boundary
    sigma = gap/2./1.89  # sigma of NGauss (best fit with SF)
    # sigma = gap/2./2.2   # sigma of NGauss (best fit with SF)
    E0    = 1.           # top of NGauss   (best fit with SF)

    z = np.linspace(zl,zr,2*anz+1)
    Ez0_tab = [(x,0.,E0*NGauss(x,sigma,0.)) for x in z]
    # display(Ez0_tab,'slice')
    poly  = Kpoly(z,sigma,0.,E0*1.)

    zstep = (zr-zl)/500.
    z = np.arange(zl,zr,zstep)
    Ez0_tab = [(x,0.,Ipoly(x, poly)) for x in z]
    display(Ez0_tab,'NG-poly')

            # OVERAP?
    # poly1 = Kpoly(z,sigma,gap,E0*1.)  # next cavity
    # Ez1_tab = [(x,0.,Ipoly(x,poly1)) for x in z]
    # display(Ez1_tab,'poly1')
    # Ez2_tab = []
    # for i in range(len(z)):
    #     x = z[i]
    #     Ez2 = Ez0_tab[i][2]+Ez1_tab[i][2]
    #     Ez2_tab.append((x,0.,Ez2))
    # display(Ez2_tab,'poly2')
            # END-OVERAP?

    # TTF calculations
    # v0  = V0(poly,zl,zr)
    tk  = Tk(poly,k,zl,zr)
    # sk  = Sk(poly,k,zl,zr)
    tkp = Tkp(poly,k,zl,zr)
    # skp = Skp(poly,k,zl,zr)
    # DEBUG('V0',v0)
    DEBUG('T(k)',tk)
    DEBUG("T'(k)",tkp)
    # DEBUG('S(k)',sk)
    # DEBUG("S'(k)",skp)

def test3(input_file):
    '''
    Second order polynomial fit with Lagrange 3 point formula to SF data
    '''
    particle = Proton(tkin=100.)
    beta     = particle.beta
    c        = PARAMS['lichtgeschwindigkeit']
    freq     = PARAMS['frequenz']
    k        = 2*pi*freq/(c*beta)

    Ez0_tab = SFdata(input_file,1.)
    # display(Ez0_tab,'slice')
    poly    = KpolySF(Ez0_tab,anz=10)
    # zl  = poly[0].zl
    # zr  = poly[-1].zr
    gap = 4.8
    zl  = -gap/2.
    zr  = +gap/2.

    zstep = (zr-zl)/500.
    z = np.arange(zl,zr,zstep)
    Ez0_tab = [(x,0.,Ipoly(x, poly)) for x in z]
    display(Ez0_tab,'SF-poly')

    # TTF calculations
    v0  = V0(poly,zl,zr)
    tk  = Tk(poly,k,zl,zr)
    # sk  = Sk(poly,k,zl,zr)
    tkp = Tkp(poly,k,zl,zr)
    # skp = Skp(poly,k,zl,zr)
    DEBUG('V0',v0)
    DEBUG('T(k)',tk)
    DEBUG("T'(k)",tkp)
    # DEBUG('S(k)',sk)
    # DEBUG("S'(k)",skp)

if __name__ == '__main__':
    input_file='SF_WDK2g44.TBL'
    ax = pre_plt(input_file)
    # test0(input_file,ax,plot=True)  # SF
    # test1()                         # NG
    test2()                         # poly-fit to NG
    test3(input_file)               # poly-fit to SF
    post_plt(ax)
