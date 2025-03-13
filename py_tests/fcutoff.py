from math import pi
from scipy.constants import c as cl

r= 0.65e-2
fc=1.8412*cl/(2*pi*r)
rc=1.8412*cl/(fc*2*pi)
# print(f'r[m] {r}, rcutoff[mm] {rc*1e3:.3e}, fcutoff[GHz] {fc*1.e-9:.3g}')
print(f'rcutoff[mm] {rc*1e3:.3e}, fcutoff[GHz] {fc*1.e-9:.3g}')