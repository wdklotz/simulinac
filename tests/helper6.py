import sys
sys.path.insert(0,'..')

from math import pi,sqrt,fabs,log
from setutil import PARAMS

c      = PARAMS['lichtgeschwindigkeit']
f0     = PARAMS['frequenz']

print('\nCutoff @ {:8.4f} [MHz] and bore R'.format(f0*1.e-6))

R     = 0.02     #bore radius
omega = 2.405 * c / R
f     = omega / 2. / pi
rhole = R

print('(R[mm],f[MHz])',(R*1.e3,f*1.e-6))

f     = f0
omega = 2.* pi * f
R     = 2.405 * c / omega
rc    = R       #cutoff

print('(R[mm],f[MHz])',(R*1.e3,f*1.e-6))

kz2 = (2.405/rc)**2 - (2.405/rhole)**2
kz  = sqrt(fabs(kz2))
z2z0 = log(100.)/kz    #daempfung 1/100 = 1%
print('propagation factor |kz|',kz)
print('penetration distance [mm] fuer E/E0=1%',z2z0*1.e3)