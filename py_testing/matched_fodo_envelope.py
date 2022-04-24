"""
    Based on Analytic Solution of the Envelope Equations for an Undepressed Matched Beam in a
    Quadrupole Doublet Channel, O.A. Anderson, L.L. LoDestro
"""
import sys
sys.path.insert(0,'..')
from setutil import Proton
from math import sqrt, sin, cos, sinh, cosh, degrees, radians, acos, pi
import numpy as np
import matplotlib.pyplot as plt

def AN(W,bgrad,z,lquad,lcav,lx,epsi,ncav):
    """
    calculates a-norm = a(z)/sqrt(epsilon*L) 
    """

    proto = Proton(tkin=W)
    brho = proto.brho
    kappa = bgrad/brho
    
    etaL = lquad
    d1 = 2*ncav*lcav+lx
    d2 = 2*ncav*lcav+lx
    L = (2*etaL+d1+d2)/2.
    d = (d1+d2)/2.
    mu = (d1-d2)/(2.*d)
    k = sqrt(kappa)
    kLpi = k*L/pi
    nu = k*d
    nu1 = nu*(1-mu)
    nu2 = nu*(1+mu)
    phi = etaL*k
    sn = sin(phi)
    cs = cos(phi)
    sh = sinh(phi)
    ch = cosh(phi)
    B = nu*ch+(1-mu**2)/2*nu**2*sh
    csigma = (ch+nu*sh)*cs - B*sn
    
    def F(phi,z):
        res= (ch+nu*sh)*sn+nu*mu*sh*sin(phi*(1.-2.*z/etaL))+B*cs+(B+sh)*cos(phi*(1-2.*z/etaL))
        return res
    try:
        a2 = F(phi,z)/sqrt(1-csigma**2)*epsi*etaL
        an = sqrt(a2)/sqrt(epsi*L)
    except ValueError:
        an = 0.
    return an

# Proton
W = 5               # kin.energy [MeV]
proto = Proton(tkin=W)
m0c2  = proto.m0c2
gb    = proto.gamma_beta
brho  = proto.brho
# Lattice
lquad = 0.04        # QUAD len
lcav  = 0.044       # CAV len
lx    = 0.1         # D between cavs
ncav  = 10          # nbof CAVs in a group
epsi  = 1.e-6       # emittance [m*rad]
# Anderson's variables
etaL = lquad
d1   = 2*ncav*lcav+lx          # space between F- and D-quad
d2   = 2*ncav*lcav+lx          # space between F- and D-quad
L    = (2*etaL + d1 + d2)/2.   # 1/2 fodo cell length
mu = (d1-d2)/2.
d  = (d1+d2)/2.
d1 = d*(1-mu)
d2 = d*(1+mu)

def Printed():
    # grads = [18,23,28,33,38,43]
    grads = [6.7]                  # F gradient [T/m]  (same for D)
    for bgrad in grads:
        kappa = bgrad/brho        # [1/m^2]
        
        # printout
        k  = sqrt(kappa)
        kLpi = k*L/pi
        nu = k*d
        nu1 = nu*(1-mu)
        nu2 = nu*(1+mu)
        phi = etaL*k
        sn = sin(phi)
        cs = cos(phi)
        sh = sinh(phi)
        ch = cosh(phi)
        B = nu*ch+(1-mu**2)/2*nu**2*sh
        csigma = (ch+nu*sh)*cs - B*sn
    
        zm  = lquad/2.              # middle of F-Quad
        an  = AN(W,bgrad,zm,lquad,lcav,lx,epsi,ncav)
        
        print("PROTON's variables\n \
            W [Mev].................: {:10.3f}\n \
            Bgrad [T/m].............: {:10.3f}\n \
            kappa [m^-2]............: {:10.3f}\n".format(W,bgrad,kappa))
        
        print("Anderson's variables\n \
            etaL....................: {:10.3f}\n \
            d1......................: {:10.3f}\n \
            d2......................: {:10.3f}\n \
            L.......................: {:10.3f}\n \
            mu......................: {:10.3f}\n \
            k.......................: {:10.3f}\n \
            kLpi....................: {:10.3f}\n \
            nu......................: {:10.3f}\n \
            nu1.....................: {:10.3f}\n \
            nu2.....................: {:10.3f}\n \
            phi.....................: {:10.3f}\n \
            B.......................: {:10.3f}\n \
            csigma..................: {:10.3f}\n".format(etaL,d1,d2,L,mu,k,kLpi,nu,nu1,nu2,phi,B,csigma))
        
        print("Results\n \
            sigma [degrees].........: {:10.3f}\n \
            a-norm..................: {:10.3f}\n".format(degrees(acos(csigma)),an))

def TwoD():
    grads = np.linspace(3,12,30)     # B-field gradients 5 to 45 T/m W=5
    # grads = np.linspace(10,20,30)    # B-field gradients 5 to 45 T/m W=25
    # grads = np.linspace(5,40,30)     # B-field gradients 5 to 45 T/m W=50
    # grads = np.linspace(5,60,30)     # B-field gradients 5 to 45 T/m W=100
    zm    = lquad/2.                 # middle of 1st quad
    abscisse = []
    ordinate = []
    for bgrad in grads:
        abscisse.append(bgrad)
        an = AN(W,bgrad,zm,lquad,lcav,lx,epsi,ncav)
        ordinate.append(an)
    
    fig,ax = plt.subplots()
    ax.set_xlabel(r"gradient    $B'[T/m]$")
    ax.set_ylabel(r"envelope    $X/\sqrt{\epsilon*L}$")
    ax.set_title("p($T_{{kin}}= {} MeV$)".format(W))
    plt.plot(abscisse,ordinate)
    # plt.show()

def ThreeD():
    # This import registers the 3D projection, but is otherwise unused.
    from mpl_toolkits.mplot3d import Axes3D

    grads = np.linspace(3,18,30)     # B-field gradients
    Wkin  = np.linspace(5,20*10,20)  # kin. energies
    zm    = lquad/2.                # FQ center
    
    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_zlim(0.7,1.5)
    ax.set_xlabel('[T/m]')
    ax.set_ylabel('[MeV]')
    ax.set_zlabel(r"$X/\sqrt{\epsilon*L}$")

    # Make data.
    X = grads
    Y = Wkin
    nx = len(X)
    ny = len(Y)
    for j in range(ny):
        z  = []
        yj = np.array([Y[j]])
        for i in range(nx):
            aij = AN(Y[j],X[i],zm,lquad,lcav,lx,epsi,ncav)
            z.append(aij)
        Z = np.array(z)
        ax.scatter(X,yj,Z)

if __name__ == '__main__':
    Printed()
    # TwoD()
    ThreeD()
    plt.show()

    
