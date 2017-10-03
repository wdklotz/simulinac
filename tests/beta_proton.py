import sys
sys.path.insert(0,'..')

from math import sqrt
from setutil import tblprnt,Proton,PARAMS

c = PARAMS['lichtgeschwindigkeit']
def BetaGamma(particle):
    vclassic   = sqrt(2.*particle.tkin/particle.e0*c**2)
    beta       = particle.beta
    gamma      = particle.gamma
    gamma_beta = particle.gamma_beta
    pc         = gamma_beta*particle.e0
    return (vclassic, beta, gamma,gamma_beta,pc)

print ('Proton Ruhemasse {:.3f} [Mev/c**2]'.format(Proton().e0))
header = [
    'T[MeV]',
    'v(class)[m/s]',
    'beta',
    'v(rel)[m/s]',
    'gamma',
    'gam*bet',
    'pc[MeV]'
    ]
rows   = []
for tk in range(5,21,1):
    vclassic,beta,gamma,gamma_beta,pc = BetaGamma(Proton(tkin=tk))
    row = [
        '{:.1f}'.format(tk),
        '{:10.4e}'.format(vclassic),
        '{:.4f}'.format(beta),
        '{:10.4e}'.format(beta*c),
        '{:1.4f}'.format(gamma),
        '{:1.4f}'.format(gamma_beta),
        '{:1.4f}'.format(pc)
        ]
    rows.append(row)
print(tblprnt(header,rows))
rows   = []
for tk in range(20,205,5):
    vclassic,beta,gamma,gamma_beta,pc = BetaGamma(Proton(tkin=tk))
    row = [
        '{:.1f}'.format(tk),
        '{:10.4e}'.format(vclassic),
        '{:.4f}'.format(beta),
        '{:10.4e}'.format(beta*c),
        '{:1.4f}'.format(gamma),
        '{:1.4f}'.format(gamma_beta),
        '{:1.4f}'.format(pc)
        ]
    rows.append(row)
print(tblprnt(header,rows))
