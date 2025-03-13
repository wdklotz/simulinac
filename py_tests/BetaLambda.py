import scipy.constants as K
from math import sqrt


def main(T):
    c  = K.c # [m/s]
    q  = K.e # [coulomb]
    mp = K.value('proton mass energy equivalent in MeV')
    me = K.value('electron mass energy equivalent in MeV')
    freq = 750e+6 # 1/sec

    E0    = mp
    beta  = sqrt(1-1/(1+T/mp)**2)
    gamma = 1+T/mp
    lamb  = c/freq
    betlamb = beta*lamb
    print(f'T {T}, beta {beta:.4f}, gamma {gamma:.4f}, freq {freq*1e-6:.0f}, lambda {lamb:.4f}, betalamda {betlamb:.4f}, 1/2*betalamda {1/2*betlamb:.4f}, 3/8*betalamda {3/8*betlamb:.4f}')

main(5)
main(10)
main(20)
main(50)
main(100)
main(200)
