import sys
sys.path.insert(0,'..')
import numpy as np
import math as M
import matplotlib.pyplot as plt
import sympy as S
from setutil import Twiss

S.init_printing()

def DEBUG_ON(*args):
    print(args)
def DEBUG_OFF(*args):
    pass
    
my_debug = DEBUG_ON
my_debug = DEBUG_OFF

def Mfodo(f,L): # (6.1)   full cell matrix
    """ 
    Wiedemann corrected
    L = distance between quad centers
    f = focus of 1/2 quad
    """
    fstar = 1./(2.*L/f**2*(1-L/f))      # Wiedemann (6.1) is wrong!!
    fstar = 1./(2.*L/f**2*(1+L/f))      # corrected
    m11 = 1.-2.*L**2/f**2
    m12 = 2.*L*(1.-L/f)
    m21 =  -1./fstar
    m22 = m11
    M = np.array([[m11,m12],[m21,m22]])
    return M

def Md(L):
    return np.array([[1.,L],[0.,1.]])        # drift
def Mq(f,L):
    return np.array([[1.,0.],[-1./f,1.]])    # thin quad
def MQ(f,L):
    k = 1./(f*L)
    if k > 0:
        phi = M.sqrt(k)*L
        m11 = m22 = M.cos(phi)
        m12 = M.sin(phi)/M.sqrt(k)
        m21 = -M.sqrt(k)*M.sin(phi)
    elif k < 0.:
        phi = M.sqrt(abs(k))*L
        m11 = m22 = M.cosh(phi)
        m12 = M.sinh(phi)/M.sqrt(abs(k))
        m21 = M.sqrt(abs(k))*M.sinh(phi)
    return np.array([[m11,m12],[m21,m22]])
def mmult(m1,m2):
    m1m2 = np.array([[0.,0.],[0.,0.]])
    for i in range(2):
        for k in range(2):
            for j in range(2):
                m1m2[i,k] += (m1[i,j] * m2[j,k])
    return m1m2
def twmatrix(M):
    C  = M[0][0]
    S  = M[0][1]
    CP = M[1][0]
    SP = M[1][1]
    m11 = C**2
    m12 = -2.*S*C
    m13 = S**2
    m21 = -C*CP
    m22 = SP*C+S*CP
    m23 = -S*SP
    m31 = CP**2
    m32 = -2.*SP*CP
    m33 = SP**2
    twissmatrix = np.array([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])
    return twissmatrix

# Bg=30.            # [T/m] DB/Dz
Bg=21.575         # [T/m] DB/Dz
E0=0.93826        # [GeV] proton rest
Ql=0.04           # full quad [m]
R=0.011           # aperture [m]

def test0():
    """ Wiedemann's FODO Formeln """
    print('--------------------------------------------" Wiedemann`s FODO Formeln "')
    print(F'kappa-opt = {M.sqrt(2):8.4g}   phi-opt = {90.} [deg]')
    print('   {T:8s}    {L:8s}  {k:8s}  {beta:8s}    {betap:8s}   {betam:8s}   {kappa:8s}   {phi:8s}'.format(T='T[MeV]',L='L[m]',k='k[m^-2]',beta='beta-kin',betap='beta+', betam='beta-', kappa='kappa', phi='phi [deg]'))
    L=0.596  # 1/2 cell 
    vrs = [(x,L) for x in [25.,22.5,20.,17.5]]
    L=0.35
    vrs += ([(x,L) for x in [15.,10.,7.5,5.]])
    for T,L in vrs:
        
        T=T*1.e-3         # [GeV] kin.energy 
        betakin=M.sqrt(1.-(1+T/E0)**-2)   # beta kinetic
        E=E0+T           # [GeV] energy
        k=0.2998*Bg/(betakin*E)    # [1/m^2]
        f=1./(k*Ql)
            
        kappa       = f/L      # (6.4) mus[50.,40.,30.,20.,15.]t be > 1
        betap       = L*kappa*(kappa+1)/M.sqrt(kappa**2-1.)   # (6.3) beta+
        betam       = L*kappa*(kappa-1)/M.sqrt(kappa**2-1.)   # (6.5) beta-
        cosphi0     = (kappa**2-2)/kappa**2   # (6.7)
        phi0        = M.acos(cosphi0)   #phi-opt
        phi0_deg    = M.degrees(phi0)
        kappa_opt   = M.sqrt(2.)         # kappa-opt for round beam   (6.11)
        phi_opt_deg = 90.                # phi-opt for round beam   (6.12)
        betap2beta_opt = kappa*(kappa+1.)/((2.+M.sqrt(2.))*M.sqrt(kappa**2-1.))    # (6.16) beta/beta-opt
        betam2beta_opt = kappa*(kappa-1.)/((2.-M.sqrt(2.))*M.sqrt(kappa**2-1.))    # (6.16)
        
        Dkappa    = abs(kappa - kappa_opt)       # small for best
        Dphi_deg = abs(phi_opt_deg - phi0_deg)   # small for best
        epsimx = R**2/(4.*L)                     # (6.15) max emittance
        
        print('{T:8.4g}   {L:8.4g}   {k:8.4g}   {beta:8.4g}   {betap:8.4g}   {betam:8.4g}   {kappa:8.4g}   {phi:8.4g}'.format(T=T*1.e3, L=L, k=k, beta=betakin, betap=betap, betam=betam, kappa=kappa, phi=phi0_deg))
    return

def test1():
    """
    Wiedemann's Formel hat einen Fehler!
    """
    print("\n\n...................",    """Wiedemann's Formel hat einen Fehler!""")
    fs,ld = S.symbols('fs ld')
    QF  = S.Matrix([[1,0],[-1/fs,1]])
    DR  = S.Matrix([[1,ld],[0,1]])
    QD  = S.Matrix([[1,0],[1/fs,1]])
    FODO1 = QD*DR*QF
    FODO2 = FODO1.subs(fs, -fs)
    
    FODO = FODO1*FODO2
    FODO = S.expand(FODO)
    print('... this is the correct one:')
    S.pprint(FODO)

    T = 25.           # kin. energy [MeV]
    T=T*1.e-3         # [GeV] kin.energy 
    betakin=M.sqrt(1.-(1+T/E0)**-2)   # beta kinetic
    E=E0+T            # [GeV] energy
    k=0.2998*Bg/(betakin*E)    # [1/m^2]
    L=0.596                    # distance between quads [m]
    Ql=0.04                    # full quad length [m]
    f = 1./(k*Ql)
    fodo1 = FODO.subs(fs,2*f).subs(ld,L/2.)
    fodo1 = np.array(fodo1)
    fodo2 = Mfodo(2*f,L/2.)         # Wiedemann  (the book is wrong!!)
    my_debug('matrix probe: fodo1-fodo2 must be zero matrix')
    zero = fodo1-fodo2
    my_debug(f'{abs(zero[0][0])}  {abs(zero[0][1])}')
    my_debug(f'{abs(zero[1][0])}  {abs(zero[1][1])}')
    return
    
def test2():
    """
    To understand and to confirm the correctnes of the SIMULINAC plots,
    check the influence of twiss-alha parameter using a thin lens FODO.
    """
    print('test2 ............................................................',
    """
    To understand and to confirm the correctnes of the SIMULINAC plots,
    this test checks the influence of the twiss-alpha parameter using a 
    thin lens FODO.
    """)
    # initial conditions for twiss parameters
    beta = 1.          # [m]
    alfa = -0.5
    epsi = 1.e-6       # emittance [m*rad]
    # other initials
    NC=56             # NC cells
    T = 25.           # kin. energy [MeV]
    T=T*1.e-3         # [GeV] kin.energy 
    betakin=M.sqrt(1.-(1+T/E0)**-2)   # beta kinetic
    E=E0+T            # [GeV] energy
    k=0.2998*Bg/(betakin*E)    # [1/m^2]
    L=0.596                    # distance between quads [m]
    Ql=0.04                    # full quad length [m]
    f=1./(k*Ql)                # ff = -fd = f
    mfodo1 = Mq(-2*f,Ql/2.)    # 1/2 qd
    mfodo2 = mmult(mfodo1,Md(L/2.))
    mfodo3 = mmult(mfodo2,Mq(2*f,Ql/2.))   # 1/2 qf
    m11 = mfodo3[0][0]
    m12 = mfodo3[0][1]
    m21 = mfodo3[1][0]
    m22 = mfodo3[1][1]
    mfodor = np.array([[m22,m12],[m21,m11]])  # reverse
    
    mfodo = mmult(mfodo3,mfodor)
    trace = abs(mfodo.trace())
    
    mfodo2 = Mfodo(2*f,L/2.)         # Wiedemann
    my_debug('matrix probe: mfodo-mfodo2 must be zero matrix')
    zero = mfodo-mfodo2
    my_debug(f'{abs(zero[0][0]):.5f}  {abs(zero[0][1]):.5f}')
    my_debug(f'{abs(zero[1][0]):.5f}  {abs(zero[1][1]):.5f}')
    
    # now in slices
    slices = 20
    ld = L/slices
    lq = Ql/slices
    f = 1./(k*lq)
    md   = Md(ld/2)         # element
    mqf  = Mq(+2*f,lq/2)    # element
    mqd  = Mq(-2*f,lq/2)    # element
    twd  = twmatrix(md)     # twiss-matrix
    twqf = twmatrix(mqf)    # twiss-matrix
    twqd = twmatrix(mqd)    # twiss-matrix
    lattice   = []          # element list
    twlattice = []          # twiss matrix list
    pos = [0.]              # position list
    s = 0.                  # postion (abzsisse)
    mfodo5 = np.array([[1.,0.],[0.,1.]])# unit matrix
    for m in range(NC):
        mfodo6 = mfodo5.dot(mfodo2)     # chain Wiedemann cells
        mfodo5 = mfodo6
        for n in range(1,slices+1):     # chain slices
            s += lq
            pos.append(s)
            lattice.append(mqd)
            twlattice.append(twqd)
        for n in range(1,slices+1):
            s += ld
            pos.append(s)
            lattice.append(md)
            twlattice.append(twd)
        for n in range(1,slices+1):
            s += lq
            pos.append(s)
            lattice.append(mqf)
            twlattice.append(twqf)
        for n in range(1,slices+1):
            s += lq
            pos.append(s)
            lattice.append(mqf)
            twlattice.append(twqf)
        for n in range(1,slices+1):
            s += ld
            pos.append(s)
            lattice.append(md)
            twlattice.append(twd)
        for n in range(1,slices+1):
            s += lq
            pos.append(s)
            lattice.append(mqd)
            twlattice.append(twqd)
    # full lattice matrix
    mfodo3 = np.array([[1.,0.],[0.,1.]])
    for node in lattice:
        mfodo4 = mfodo3.dot(node)
        mfodo3 = mfodo4
    my_debug('lattice probe: mfodo3-mfodo5 must be zero matrix')
    zero = mfodo3-mfodo5
    my_debug(f'{abs(zero[0][0]):.5f}  {abs(zero[0][1]):.5f}')
    my_debug(f'{abs(zero[1][0]):.5f}  {abs(zero[1][1]):.5f}')
    
    twx = Twiss(beta,alfa,epsi)
    # C track
    track_point_C = np.array(twx.y1())
    # S track
    track_point_S = np.array(twx.y4())

    points_c      = [track_point_C]
    points_s      = [track_point_S]
    for node in lattice:
        C_point = node.dot(track_point_C)
        points_c.append(C_point)
        track_point_C = C_point
        S_point = node.dot(track_point_S)
        points_s.append(S_point)
        track_point_S = S_point

    # enveloppe 
    beta0, alfa0, gamma0 , epsi0 = twx()
    twiss_point = np.array([beta0,alfa0,gamma0])   # twiss track
    twpoints = [twiss_point]
    for bmx in twlattice:
        twissp = bmx.dot(twiss_point)
        twpoints.append(twissp)
        twiss_point = twissp

    # trajectories
    c_trk  = [v[0] for v in points_c]
    s_trk  = [v[0] for v in points_s]
    # enveloppe  (beta*emittance)^(1/2)
    sgx = [M.sqrt(v[0]*epsi0) for v in twpoints]
    
    plt.figure('x & sgx')
    ax1 = plt.subplot(111)
    ax1.set_title(F'alpha = {alfa:2.1f}')
    ax1.plot(pos,c_trk,   label='C')
    ax1.plot(pos,s_trk,   label='S')
    ax1.plot(pos,sgx,     label='env')
    ax1.legend(loc='lower right',fontsize='small')
    plt.show()
    return
    
if __name__ == '__main__':
    test0()
    test1()
    test2()