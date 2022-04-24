import sys
sys.path.insert(0,'..')

from math import sqrt,pi
from setutil import tblprnt,dictprnt
import scipy.constants as C

# Vorgaben
# freq [Hz], T [Mev}, emitw [rad], DT2T=delta-T/T [%]
freq  = 816.0e6
T     = 50
emitw = 1.34e-5
DT2T  = 2.e-1

# ------------------------------------
DT2T   = DT2T*1.e-2    # [%] -> [1]
clight = C.c
m0c2   = C.value('proton mass energy equivalent in MeV')
lamb   = clight/freq       # wavelength
gamma  = 1.+T/m0c2         # gamma lorentz
gabe   = sqrt(gamma**2-1.) # gamma*beta lorentz
beta   = gabe/gamma        # beta lorentz
bela   = beta*lamb         # beta*lambda
Dgamma = (gamma-1.)*DT2T   # delta-gamma==w
sigphi = emitw/Dgamma
z      =-bela/(2.*pi)*sigphi # [m]
Dp2p   = 1./(gabe*beta)*Dgamma

dicto = {}
# * means Vorgabe
dicto['frequency [Hz]*']               = '{:8.3e} freq'.format(freq)
dicto['kinetic energy T [MeV]*']       = '{:8.3e} T'.format(T)
dicto['emit{phi-w} [rad]*']            = '{:8.3e} emitw'.format(emitw)
dicto['delta-T/T spread*']             = '{:8.3e} DT2T'.format(DT2T)
dicto['gamma Lorentz']                 = '{:8.3e} gamma'.format(gamma)
dicto['beta Lorentz']                  = '{:8.3e} beta'.format(beta)
dicto['wavelength [m]']                = '{:8.3e} lamb'.format(lamb)
dicto['phi spread [rad]']              = '{:8.3e} Dphi'.format(sigphi)
dicto['z spread [m]']                  = '{:8.3e} bunch'.format(z)
dicto['w spread']                      = '{:8.3e} Dgamma, dE/E0'.format(Dgamma)
dicto['Dp/p spread']                   = '{:8.3e} Dp2p'.format(Dp2p)


print(dictprnt(dicto,"parameter for BEAM Input"))











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