import sys
sys.path.append('..')
# print(sys.path)

from math import sqrt
from setutil import tblprnt,Proton,PARAMS
import scipy.constants as C

c = C.c             #[m/s]
freq = 816. if len(sys.argv) == 1 else float(sys.argv[1])
freq = freq*1.e6    #[Hz]

def BetaGamma(particle):
    vclassic   = sqrt(2.*particle.tkin/particle.e0*c**2)
    beta       = particle.beta
    gamma      = particle.gamma
    gamma_beta = particle.gamma_beta
    pc         = gamma_beta*particle.e0
    brho       = gamma_beta*3.1297
    lamb       = c/freq
    blamb      = beta*lamb*100.
    return (vclassic, beta, gamma, gamma_beta, pc, brho, blamb)

print ('Proton rest mass {:.3f} [Mev/c**2],  frequency[MHz] {:.3f}'.format(Proton(50.).e0,freq*1.e-6))
header = [
    'T[MeV]',
    'v(class)[m/s]',
    'beta',
    'v(rel)[m/s]',
    'gamma',
    'gam*bet',
    'pc[MeV]',
    'B*rho[Tm]',
    'beta*lambda[cm]'
    ]
rows   = []
for tk in range(5,21,1):
    vclassic, beta, gamma, gamma_beta, pc, brho, blamb = BetaGamma(Proton(tkin=tk))
    row = [
        '{:.1f}'.format(tk),
        '{:10.4e}'.format(vclassic),
        '{:.4f}'.format(beta),
        '{:10.4e}'.format(beta*c),
        '{:1.4f}'.format(gamma),
        '{:1.4f}'.format(gamma_beta),
        '{:1.4f}'.format(pc),
        '{:1.4f}'.format(brho),
        '{:1.4f}'.format(blamb)
        ]
    rows.append(row)
print(tblprnt(header,rows))
rows   = []
for tk in range(20,205,5):
    vclassic, beta, gamma, gamma_beta, pc, brho, blamb = BetaGamma(Proton(tkin=tk))
    row = [
        '{:.1f}'.format(tk),
        '{:10.4e}'.format(vclassic),
        '{:.4f}'.format(beta),
        '{:10.4e}'.format(beta*c),
        '{:1.4f}'.format(gamma),
        '{:1.4f}'.format(gamma_beta),
        '{:1.4f}'.format(pc),
        '{:1.4f}'.format(brho),
        '{:1.4f}'.format(blamb)
        ]
    rows.append(row)
print(tblprnt(header,rows))
