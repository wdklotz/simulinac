import numpy as np
import math as M

def Mfodo(f,L): # (6.1)   cell matrix
    fstar = 1./(2*(1.-L/f)*(L/f**2))
    m11 = 1.-2.*L**2/f**2
    m12 = 2.*L*(1.+L/f)
    m21 =  -1./fstar
    m22 = m11
    M = np.array([[m11,m12],[m21,m22]])
    return M

# g=30.            # [T/m] DB/Dz
g=21.575         # [T/m] DB/Dz
E0=0.93826       # [GeV] proton rest

print('   {0:8s}    {1:8s}    {2:8s}  {3:8s}   {4:8s}  {5:8s}  {6:8s}   {7:8s}'.format('T[MeV]','L[m]','k[m^-1]','beta rel.','beta+', 'beta-', 'kappa', 'phi [deg]'))
L=0.596# 1/2 cell 
vrs = [(x,L) for x in [25.,22.5,20.,17.5]]
L=0.35
vrs += ([(x,L) for x in [15.,10.,7.5,5.]])
for T,L in vrs:
    Ql=0.04          # quad [m]
    R=0.011          # aperture [m]
    # d3l=0.080        # d3 [m]
    # d5l=0.022        # d5 [m]
    # N=9              # N gaps
    
    T=T*1.e-3         # [GeV] kin.energy 
    betein=M.sqrt(1.-(1+T/E0)**-2)   # beta relativistic
    E=E0+T           # [GeV] energy
    k=0.2998*g/(betein*E)    # [1/m^2]
    f=1./(k*Ql)
        
    mf = Mfodo(f,L)
    md = Mfodo(-f,L)
    
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
    
    print('{0:8.4g}   {1:8.4g}   {2:8.4g}   {3:8.4g}   {4:8.4g}   {5:8.4g}   {6:8.4g}   {7:8.4g}'.format(T*1.e3, L, k, betein, betap, betam, kappa, phi0_deg))
    pass