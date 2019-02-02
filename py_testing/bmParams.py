import sys
sys.path.insert(0,'..')

from math import sqrt,pi
from setutil import tblprnt,dictprnt
import scipy.constants as C

# Vorgaben
# freq [Hz], W [Mev}, emitw [rad], DW2W=delta-W/W []
freq  = 816.0e6
W     = 25
emitw = 1.e-6
DW2W  = 1.e-2

# ------------------------------------
clight = C.c
m0c2   = C.value('proton mass energy equivalent in MeV')
lamb   = clight/freq       # wavelength
gamma  = 1.+W/m0c2         # gamma lorentz
gabe   = sqrt(gamma**2-1.) # gamma*beta lorentz
beta   = gabe/gamma        # beta lorentz
bela   = beta*lamb         # beta*lambda
Dgamma = (gamma-1.)*DW2W   # delta-gamma==w
sigphi = emitw/Dgamma
z      =-bela/(2.*pi)*sigphi # [m]

dicto = {}
# * means Vorgabe
dicto['frequency: f[Hz]*']                 = freq
dicto['kinetic energy: W[MeV]*']           = W
dicto['emittance{phi-w}: emitw[rad]*']     = emitw
dicto['kin.energy spread=deltaW/W: DW2W*'] = DW2W
dicto['gamma Lorentz']                     = gamma
dicto['beta Lorentz']                      = beta
dicto['wavelength: lamb[m]']               = lamb
dicto['sigma_phase: sigphi[rad]']          = sigphi
dicto['sigma_z: z[m]']                     = z
dicto['sigma_e=w=deltaW/m0c2: Dgamma']     = Dgamma


print(dictprnt(dicto,"parameter for BMAD"))











# def delta_phi(T,z):
#     gamma = (T+m0)/m0
#     beta  = sqrt(1.-1./gamma**2)
#     phi   = 360/beta/lamb*z
#     return phi
# 
# header = ['T[Mev]','z[m]','delta phi[deg]','delta phi[deg/cm]']
# tablrows = []
# for t in [10.,50.,100.,150.,200.,250.,1.e3]:
#     dphi = delta_phi(t,z)
#     tablrow = ['{:8.4g}'.format(t),
#                 '{:8.4f}'.format(z),
#                 '{:8.4f}'.format(dphi),
#                 '{:8.4f}'.format(dphi/(z*1.e2))]
#     tablrows.append(tablrow)
# print('Phase advance as function of kinetic energy T and distance z')
# print('c={:8.4g}[m/s], m0={:8.4g}[MeV/c2], fr={:8.4g}[MHz]'.format(c,m0,freq*1.e-6))
# print(tblprnt(header,tablrows))